#!/bin/bash
#SBATCH -J _SJRN_
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=_nCPUs_
#SBATCH --qos=gp_debug
#SBATCH --time=0:10:00
#SBATCH --account=bsc72

ml intel/2023.2.0
ml cmake/3.25.1
ml impi/2021.10.0
ml mkl/2023.2.0
ml miniconda/24.1.2

source activate /gpfs/projects/bsc72/conda_envs/nbdsuite/0.5.0
umask u=rwx,g=rwx,o=

python -m nbdsuite.main nbdsuite_input.yaml
