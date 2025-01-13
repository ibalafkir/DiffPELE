#!/bin/bash

# Description: From a PELE equilibration simulation run using several initial systems
# (40), get the pose with best total energy (TE) among the top 20% in binding energy
# (BE)

# Usage: bash script.sh /path/to/PELE/dir

# -----------------------------------------------------------------------------------

set -e # Terminate if error

# Verifying if the input is a directory
if [ $# -ne 1 ]; then
    echo "Use: $0 directory"
    exit 1
fi

# Function to merge reports
merge() {
    # Creating a new directory
    PELE_DIR="$1"
    OUT="$PELE_DIR/merged_reports"
    mkdir -p "$OUT"

    # Merging reports
    for ((i = 0; i < 40; i += 1))
    do
        initial_system="$OUT/system_$((i+1)).txt"
        cat "$PELE_DIR/outEQ/0/report_$((i + 1))" \
            "$PELE_DIR/outEQ/0/report_$((i + 41))" \
            "$PELE_DIR/outEQ/0/report_$((i + 81))" \
            "$PELE_DIR/outEQ/0/report_$((i + 121))" \
            "$PELE_DIR/outEQ/0/report_$((i + 161))" > "$initial_system"
    done
}

# Function to extract a specific model from a trajectory PDB file
extract_model() {
    local TRAJECTORY_PDB="$1"
    local N_MODEL="$2"
    local PELE_DIR="$3"

    local TRAJECTORY_PDB_name=$(basename "${TRAJECTORY_PDB%.*}")
    local OUT="${PELE_DIR}/pdbsEQ/${TRAJECTORY_PDB_name}_model_${N_MODEL}.pdb"

    awk -v num="$N_MODEL" '
    $1 == "MODEL" && $2 == num { flag=1; next }
    flag && $1 == "ENDMDL" { flag=0; next }
    flag { print }
' "$TRAJECTORY_PDB" > "$OUT"
}

# Parameters
PELE_DIR=$1

# Merge reports
merge $PELE_DIR

# Iterate over grouped reports, get a equilibrated snapshot per grouped reports/trajectory
# and write metrics in
for i in {1..40}
do
awk '$1 == 1' "$PELE_DIR/merged_reports/system_${i}.txt" > snapshots.txt # remove headers
no_lines=$(cat snapshots.txt | wc -l)
#echo $no_lines
no_lines_filtrated=$(echo "$no_lines * 0.2" | bc)
no_lines_filtrated_aprox=$(echo "$no_lines * 0.2" | bc | awk '{print int($1+1)}') # picks up 20%+1pose
#echo $no_lines_filtrated
#echo $no_lines_filtrated_aprox
sort -n -k 5,5 snapshots.txt | head -n $no_lines_filtrated_aprox > snapshots_f.txt # get top 20% BE 
sort -n -k 4,4 snapshots_f.txt | head -n 1 > snapshotOfInterest.txt # order by TE and get top 1 TE
cat snapshotOfInterest.txt >> equilibratedSnapshots.txt # write snapshot
rm snapshotOfInterest.txt
rm snapshots_f.txt
rm snapshots.txt
done

# Debug (verify there are 40 snapshots in the equilibration file)
if [ $(wc -l < equilibratedSnapshots.txt) -ne 40 ]; then
    echo "Expected 40 equilibrated snapshots and found a different number"
    break
fi

# Filtering the equilibrated snapshots
# currently: only keeping top 33% in TE of the 40 equilibrated systems (i.e. 13)
sort -n -k 4,4 equilibratedSnapshots.txt | head -n 13 > equilibratedSnapshots_F.txt

# Getting equilibrated PDB models by looking for the snapshot using values from the last 3 metrics (V6-8)

OUT="equilibratedSnapshots_F.txt"

while IFS= read -r line
do
  acceptedStep=$(echo $line | awk '{print $3}')
  N_MODEL=$((acceptedStep + 1))
  V6=$(echo "$line" | awk '{print $6}')
  V7=$(echo "$line" | awk '{print $7}')
  V8=$(echo "$line" | awk '{print $8}')
  grep -rl "${V6}    ${V7}    ${V8}" "${PELE_DIR}/outEQ/0" --include="rep*" > snapshotOfInterest_path.txt
  
  # Debug (verify there is one only snapshot)
  if [ $(wc -l < snapshotOfInterest_path.txt) -gt 1 ]; then
	  echo "More than 1 snapshot coincided with the pattern"
	  break
  fi

  # Get the trajectory number
  N_TRAJ=$(sed 's/.*_//' snapshotOfInterest_path.txt)
  
  # Calling the function to extract the model
  extract_model "${PELE_DIR}/outEQ/0/trajectory_${N_TRAJ}.pdb" $N_MODEL $PELE_DIR
  rm snapshotOfInterest_path.txt
done < "$OUT"


rm equilibratedSnapshots.txt equilibratedSnapshots_F.txt
rm -r $PELE_DIR/merged_reports
