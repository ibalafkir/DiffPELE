# DiffPELE
This pipeline optimizes flexible protein-protein and antibody-antigen interactions. This is done through interface backbone diffusion and Monte Carlo rotations/translations using RFdiffusion (two fist gifs) and PELE (last picture). 

The objective is to accurately model the conformations of the regions involved in the antibody-antigen or protein-protein interactions.

## Installation
### Environment manual setup
```bash
git clone https://github.com/ibalafkir/DiffPELE.git
conda create --name diffpele python=3.11
conda activate diffpele
```
For developers:
```bash
python setup.py develop
```
For users:
```bash
python setup.py install
```
Install dependencies:
```bash
conda install biopandas -c conda-forge
pip install pdb-tools
pip install pyrosetta-installer
python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'
pip install DockQ
```

## Usage
The objective of this repository is to generate RFdiffusion, FastRelax and PELE runners for our cluster MareNostrumV. The idea is to execute the scripts with the required inputs to generate these runners and then run them following a particular order so as to try to optimize the interaction interface between the two protein components. Soon a guide explaining the order and steps will be indicated.

### System preparation

Given your protein-protein system, keep only amino acids from chains of interest and remove waters, small molecules, etc.
Use the script ./DiffPELE/src/0_preprocess_system.py to fastly do PDB operations (remove/save chains, remove hetatms, add TERs, renumber antibody insertion numbers...). The prepared PDB must not have residue numbers gaps or antibody insertion codes (1, 1A, 1B, 2...).

Finally, run this Schrodinger pipeline:
```bash
$SCHRODINGER/utilities/prepwizard input.pdb output.pdb -rehtreat -disulfides -fillloops -fillsidechains -propka_pH 7.4 -minimize_adj_h -f OPLS_2005
```
### Energy baseline
The energy baseline corresponds to energy calculations (total and binding energy) of the conformations obtained from the conformational space that occur in the given binding mode. A protein-protein complex in which interactions are not optimized, it is expected that the ensemble of these conformations show high energy values. 

Here, we do the baseline using PELE in a protocol that consists of two stages: equilibration and production. You can compare the PELE energy profiles of an optimized and and a not-optimized version of a system.

To generate the PELE MNV runners:
```bash
python 3_pele_simulation.py -pdb ./DiffPELE/examples/4POU_b.pdb -rc B -lc A -dc 12.0 -m single
```
Run this sequentally in MNV:
```bash
sbatch nbdsuite.sh # Prepare system for PELE
sbatch runEq.sh # Run equilibration
./DiffPELE/bin/equilibrate_single.sh /path/to/PELE/directory # Equilibrate pose (top total energy in top 20% binding energy)
sbatch runProd.sh # Run production
```
Plot results in MNV modifying and using this script, that uses [AdaptivePELE](https://github.com/BSC-CNS-EAPM/AdaptivePELE)
```bash
#!/bin/bash
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
	-o /path/to/PELE/directory/outPROD/0 \
	-t scatter \
	-z _ -y 5 -x 4 --hide_logo \
	-s /path/to/save/png \
	--title '' \
	--xlowest _ --xhighest _ --ylowest _ --yhighest _ --zlowest _ --zhighest _
```

### Structural baseline
Run DockQ if you are comparing an optimized system and a non-optimized system.

```bash
DockQ <a> <b>
```
where 'a' stands for the system with non-optimized interactions and 'b' stands for the optimized one.

### The DiffPELE pipeline
Given a protein-protein system in the approx. correct binding mode, the idea is to explore the interface backbone conformations using RFdiffusion, minimize the pose while repacking side chains and adjust the binding mode and side chain orientations using PELE. This is 1 epoch: we advice running several (the more epochs are run, the more the interface backbone conformational space is explored).

#### Interface diffusion
Generate runners for RFdiffusion:
```bash
python 1_interface_diffusion.py -m setup -pdb ./DiffPELE/examples/4POU_b.pdb -rc B -lc A -dc 12.0 -rip /gpfs/projects/bsc72/Repos/RFdiffusion/scripts/run_inference.py
```
Now run the following to correct the PDB format of diffusion models:
```bash
python 1_interface_diffusion.py -m fix -pdb ./DiffPELE/examples/4POU_b.pdb -dmi /path/to/*diffModels
```

#### Minimization and side chains repacking
Generate FastRelax runner on the PDBs located in the diffModels_fix folder:
```bash
python 2_fastrelax.py -m parallel -pdb /path/to/pdb
```
and launch runner:
```bash
sbatch fastrelax.sh
```

#### Binding mode adjustment
Finally, run PELE for all the fixed and minimized diffusion models. To create the runners: 
```bash
python 3_pele_simulation.py -pdb ./DiffPELE/examples/4POU_b.pdb -rc B -lc A -dc 12.0 -m multiple -mpdb /path/to/dir/with/diffusion_models
```
And run:
```bash
sbatch nbdsuite.sh # Prepare system for PELE
sbatch runEq.sh # Run equilibration
./bin/equilibrate_multiple.sh /path/to/PELE/directory # Equilibrate pose (top total energy in top 20% binding energy)
sbatch runProd.sh # Run production
```
### Analysis
Run DockQ and compare PELE energy profiles as mentioned before.