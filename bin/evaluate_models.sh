#!/bin/bash

# Description: From a PELE production simulation, get the poses with 
# low binding energies (BE) and compare to a native system using DockQ

# Usage: bash script.sh /path/to/PELE/dir /path/to/native.pdb receptor_chains ligand_chains

# ---------------------------------------------------------------

set -e # Terminate if error


# Function to extract a specific model from a trajectory PDB file
extract_model() {
    local TRAJECTORY_PDB="$1"
    local N_MODEL="$2"
    local PELE_DIR="$3"

    local TRAJECTORY_PDB_name=$(basename "${TRAJECTORY_PDB%.*}")
    local OUT="${PELE_DIR}/structure_evaluation/${TRAJECTORY_PDB_name}_model_${N_MODEL}.pdb"

    awk -v num="$N_MODEL" '
    $1 == "MODEL" && $2 == num { flag=1; next }
    flag && $1 == "ENDMDL" { flag=0; next }
    flag { print }
' "$TRAJECTORY_PDB" > "$OUT"
}


# Parameters
PELE_DIR=$(realpath "$1")
NATIVE_PDB=$(realpath "$2")
RECEPTOR_CHAINS=$3
LIGAND_CHAINS=$4

# Make output dir
mkdir -p $PELE_DIR/structure_evaluation

# Group directories and make output file
cat $PELE_DIR/outPROD/0/rep* | awk '$5<-40' | sort -k5,5n | uniq | head -n 300
cat $PELE_DIR/outPROD/0/rep* | awk '$1 == 1' | sort -k5,5n | uniq | head -n 300 > $PELE_DIR/structure_evaluation/reports.txt
echo "You will be evaluating this number of poses:"
cat $PELE_DIR/structure_evaluation/reports.txt | wc -l
NUM_MODEL=1

touch $PELE_DIR/structure_evaluation/dockq.txt

while IFS= read -r accepted_step
do

# Debug
echo "Computing model number: ${NUM_MODEL}"

# Use accepted_step's metrics to get model to get
acceptedStep=$(echo $accepted_step | awk '{print $3}')
n_model=$((acceptedStep + 1))

# Use accepted_step's metrics to find model to get path
V6=$(echo "$accepted_step" | awk '{print $6}')
V7=$(echo "$accepted_step" | awk '{print $7}')
V8=$(echo "$accepted_step" | awk '{print $8}')
accepted_step_path=$(grep -rl "${V6}    ${V7}    ${V8}" "${PELE_DIR}/outPROD/0" --include="rep*" | head -n 1)

# Get the trajectory number
n_traj="${accepted_step_path##*_}"

# Extract model
extract_model "${PELE_DIR}/outPROD/0/trajectory_$n_traj.pdb" $n_model $PELE_DIR

# Indicate in dockq.txt that this is a new pose
echo "New pose:" >> $PELE_DIR/structure_evaluation/dockq.txt

# Write PDB path in dockq.txt
echo "${PELE_DIR}/structure_evaluation/trajectory_${n_traj}_model_${n_model}.pdb" >> $PELE_DIR/structure_evaluation/dockq.txt

# Write accepted_step's metrics in $PELE_DIR/structure_evaluation/dockq.txt
echo "$accepted_step" >> $PELE_DIR/structure_evaluation/dockq.txt

# Run DockQ
#dockq_result=$(DockQ "${PELE_DIR}/structure_evaluation/trajectory_${n_traj}_model_${n_model}.pdb" "$NATIVE_PDB" --short --mapping AB:AB | tail -n 1 | awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}') 
dockq_result=$(python /gpfs/scratch/bsc72/ismael/repos/DiffPELE/src/4_evaluate_models.py -pdb "$NATIVE_PDB" -mp "${PELE_DIR}/structure_evaluation/trajectory_${n_traj}_model_${n_model}.pdb" -rc $RECEPTOR_CHAINS -lc $LIGAND_CHAINS)

# Write DockQ output in $PELE_DIR/structure_evaluation/dockq.txt
echo $dockq_result >> $PELE_DIR/structure_evaluation/dockq.txt

# Write space	
echo " " >> $PELE_DIR/structure_evaluation/dockq.txt	

# Increase count for next model
((NUM_MODEL++))  # Increment the counter

# Print from time to time
if (( NUM_MODEL % 50 == 0 )); then
cat $PELE_DIR/structure_evaluation/dockq.txt
fi


done < $PELE_DIR/structure_evaluation/reports.txt

# Calculate baseline and DockQ and Fnats of DiffPELE models
fnats=$(cat $PELE_DIR/structure_evaluation/dockq.txt | grep 'fnat' | awk '{print $8}' | sort -nr | awk '{printf (NR==1 ? "%s" : ",%s"), $0}')
dockq=$(cat $PELE_DIR/structure_evaluation/dockq.txt | grep 'DockQ' | awk '{print $2}' | sort -nr | awk '{printf (NR==1 ? "%s" : ",%s"), $0}')
binding=$(cat $PELE_DIR/structure_evaluation/dockq.txt | awk '{print $5}' | awk '$1<0' | awk 'NF' | sort -n | awk '{printf (NR==1 ? "%s" : ",%s"), $0}')
fnats=$(cat $PELE_DIR/structure_evaluation/dockq.txt | grep 'fnat' | awk '{print $8}' | awk '{printf (NR==1 ? "%s" : ",%s"), $0}')
total=$(cat $PELE_DIR/structure_evaluation/dockq.txt | awk '{print $4}' | awk '$1<0' | awk 'NF' | sort -n | awk '{printf (NR==1 ? "%s" : ",%s"), $0}')

# Append results
echo " " >> $PELE_DIR/structure_evaluation/dockq.txt
echo " " >> $PELE_DIR/structure_evaluation/dockq.txt
echo "Fnats:" >> $PELE_DIR/structure_evaluation/dockq.txt
echo $fnats >> $PELE_DIR/structure_evaluation/dockq.txt
echo " " >> $PELE_DIR/structure_evaluation/dockq.txt
echo "DockQ:" >> $PELE_DIR/structure_evaluation/dockq.txt
echo $dockq >> $PELE_DIR/structure_evaluation/dockq.txt
echo "Binding ene:" >> $PELE_DIR/structure_evaluation/dockq.txt
echo $binding >> $PELE_DIR/structure_evaluation/dockq.txt
echo " " >> $PELE_DIR/structure_evaluation/dockq.txt
echo "Total ene:" >> $PELE_DIR/structure_evaluation/dockq.txt
echo $total >> $PELE_DIR/structure_evaluation/dockq.txt
