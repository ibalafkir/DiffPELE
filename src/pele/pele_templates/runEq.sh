#!/bin/bash
#SBATCH --job-name _EJRN_
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=_nCPUs_
#SBATCH --time=2:00:00            
#SBATCH -D .
#SBATCH --cpus-per-task=1
#SBATCH --qos=gp_bscls
#SBATCH --account=bsc72


module purge
ml intel/2023.2.0
ml cmake/3.25.1
ml impi/2021.10.0
ml mkl/2023.2.0
ml boost/1.77.0-gcc
module load anaconda

export SCHRODINGER="/gpfs/projects/bsc72/Programs/SCHRODINGER_ACADEMIC"
export SCHRODINGER_PYTHONPATH="/gpfs/projects/bsc72/Programs/SCHRODINGER_ACADEMIC/internal/lib/python2.7/site-packages"
export PELE="/gpfs/projects/bsc72/PELE++/mnv/1.8.0/bin/"
export PELE_EXEC="/gpfs/projects/bsc72/PELE++/mnv/1.8.0/bin/PELE-1.8_mpi_intel"
export SRUN=1  # this is to avoid having to set usesrun: true in input.yaml
source activate /gpfs/projects/bsc72/conda_envs/adaptive

python -m AdaptivePELE.adaptiveSampling ctrl/adEq.conf
