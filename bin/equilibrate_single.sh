#!/bin/bash

# Description: From a PELE equilibration simulation, get the pose 
# with best total energy (TE) among the top 20% in binding energy
# (BE)

# Usage: bash script.sh /path/to/PELE/dir

# ---------------------------------------------------------------

set -e # Terminate if error

# Verifying if the input is a directory
if [ $# -ne 1 ]; then
    echo "Use: $0 directory"
    exit 1
fi

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

# Group reports
for file in "$PELE_DIR"/outEQ/0/report_*
do
    cat $file >> reports.txt
done

# Remove headers
awk '$1 == 1' reports.txt > reports_snapshots.txt

# Order by binding energy
sort -k 5,5n reports_snapshots.txt | uniq > reports_snapshots_BE.txt

# Get the best 20% in binding energy
N_SNAPSHOTS=$(wc -l < reports_snapshots_BE.txt)
N_SNAPSHOTS_FILTERED=$(echo "$N_SNAPSHOTS * 0.2" | bc | awk '{print ($1 == int($1)) ? $1 : int($1)+1}' | sed 's/\..*//') # always rounds up
head reports_snapshots_BE.txt -n $N_SNAPSHOTS_FILTERED > reports_snapshots_BE_filtered.txt

# Ordering by total energy
sort -k 4,4n reports_snapshots_BE_filtered.txt > reports_snapshots_BE_filtered_TE.txt

# Get the snapshot of interest: best total energy in 20% best binding energy
head reports_snapshots_BE_filtered_TE.txt -n 1 > snapshotOfInterest.txt
echo "This is the selected pose:"
cat snapshotOfInterest.txt

# Getting equilibrated PDB model by looking for the snapshot using values from the last 3 metrics (V6-8)
acceptedStep=$(cat snapshotOfInterest.txt | awk '{print $3}')
N_MODEL=$((acceptedStep + 1)) # Models start by 1 and acceptedSteps by 0
V6=$(cat snapshotOfInterest.txt | awk '{print $6}')
V7=$(cat snapshotOfInterest.txt | awk '{print $7}')
V8=$(cat snapshotOfInterest.txt | awk '{print $8}')
grep -rl "${V6}    ${V7}    ${V8}" "${PELE_DIR}/outEQ/0" --include="rep*" | head -n 1 > snapshotOfInterest_path.txt # even if more than 1 snapshot coincides, it picks the first
echo "The selected pose lies in:"
cat snapshotOfInterest_path.txt

# Debug (verify only one pose was selected)
if [ $(wc -l < snapshotOfInterest_path.txt) -gt 1 ]; then
    echo "More than 1 snapshot coincided with the pattern"
    #break
fi

# Get the trajectory number
N_TRAJ=$(sed 's/.*_//' snapshotOfInterest_path.txt)

# Calling the function to extract the model
extract_model "${PELE_DIR}/outEQ/0/trajectory_$N_TRAJ.pdb" $N_MODEL $PELE_DIR

# Deleting temp files
rm reports.txt reports_snapshots.txt reports_snapshots_BE.txt reports_snapshots_BE_filtered.txt reports_snapshots_BE_filtered_TE.txt snapshotOfInterest.txt snapshotOfInterest_path.txt
