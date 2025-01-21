# DiffPELE
This pipeline optimizes flexible protein-protein and antibody-antigen interactions. This is done through interface backbone diffusion and Monte Carlo rotations/translations using RFdiffusion and PELE. The objective is to accurately model the conformations of the regions involved in the antibody-antigen or protein-protein interactions.

## Installation
### Environment manual setup
```bash
$ git clone https://github.com/ibalafkir/DiffPELE.git
$ conda create --name diffpele python=3.11
$ conda activate diffpele
```
For developers:
```bash
$ python setup.py develop
```
For users:
```bash
$ python setup.py install
```
Install dependencies:
```bash
$ conda install biopandas -c conda-forge
$ pip install pdb-tools
$ pip install pyrosetta-installer
$ python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'
$ pip install DockQ
```

## Usage
The objective of this repository is to generate RFdiffusion, FastRelax and PELE runners for our cluster MareNostrumV. The idea is to execute the scripts with the required inputs to generate these runners and then run them following a particular order so as to try to optimize the interaction interface between the two protein components. Soon a guide explaining the order and steps will be indicated.

### System preparation

Given your protein-protein system, keep only amino acids from chains of interest and remove waters, small molecules, etc.
Use the script ./DiffPELE/src/0_preprocess_system.py to fastly do PDB operations (remove/save chains, remove hetatms, add TERs, renumber antibody insertion numbers...). The prepared PDB must not have residue numbers gaps or antibody insertion codes (1, 1A, 1B, 2...).

Finally, run this Schrodinger pipeline:
```bash
$ $SCHRODINGER/utilities/prepwizard input.pdb output.pdb -rehtreat -disulfides -fillloops -fillsidechains -propka_pH 7.4 -minimize_adj_h -f OPLS_2005
```
### Energy baseline
The energy baseline corresponds to energy calculations (total and binding energy) of the conformations obtained from the conformational space that occur in the given binding mode. A protein-protein complex in which interactions are not optimized, it is expected that the ensemble of these conformations show high energy values.

Here, we do the baseline using PELE in a protocol that consists of two stages: equilibration and production. 

To generate the PELE MNV runners:
```bash
python 3_pele_simulation.py -pdb ../../testdiffpele/4POU_b.pdb -rc B -lc A -dc 12.0 -m single
```
Run this sequentally in MNV:
```bash
sbatch nbdsuite.sh # Prepare system for PELE
sbatch runEq.sh # Run equilibration
./bin/equilibrate_single.sh /path/to/PELE/directory # Equilibrate pose (top total energy in top 20% binding energy)
sbatch runProd.sh # Run production
```
Plot results:
```bash
module purge
ml intel/2023.2.0
ml cmake/3.25.1
ml impi/2021.10.0
ml mkl/2023.2.0
ml boost/1.77.0-gcc
ml anaconda
ml bsc
ml transfer
eval "$(conda shell.bash hook)"
conda activate /gpfs/projects/bsc72/conda_envs/platform

python -m pele_platform.plotter \
	-o /gpfs/scratch/bsc72/ismael/projects/PFL/1PXV/1PXV_ub_pele/outputPROD/0 \
	-t scatter \
	-z   -y 5 -x 4 --hide_logo \
	-s /gpfs/scratch/bsc72/ismael/1pxv_ub_pele_prod.png \
	--title '1PXV unbound production' \
	--xlowest -11700 --xhighest -11100 --ylowest -160 --yhighest 0 --zlowest 0 --zhighest 5
```


### Structural baseline
DockQ

### The DiffPELE pipeline
Explain

Mention iterations

#### Interface diffusion
RFdiffusion

#### Minimization and side chains repacking
FastRelax

#### Binding mode adjustment
PELE

### Analysis