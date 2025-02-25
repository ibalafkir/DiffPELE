#!/bin/bash

# Description: From a PELE production simulation, get the poses with 
# low binding energies (BE) and compare to a native system using DockQ

# Usage: bash script.sh /path/to/PELE/dir /path/to/native.pdb receptor_chains ligand_chains binding_threshold

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
BINDING_THRESHOLD=$5

# Make output dir
mkdir -p $PELE_DIR/structure_evaluation

# Group directories and make output file
cat $PELE_DIR/outPROD/0/rep* | awk '$1 == 1' | sort -u | awk -v threshold="$BINDING_THRESHOLD" '$5 < threshold' > $PELE_DIR/structure_evaluation/reports.txt
touch $PELE_DIR/structure_evaluation/dockq.txt

while IFS= read -r accepted_step
do

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

done < $PELE_DIR/structure_evaluation/reports.txt