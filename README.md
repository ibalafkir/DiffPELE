# DiffPELE
This pipeline optimizes flexible protein-protein and antibody-antigen interactions. This is done through interface backbone diffusion and Monte Carlo rotations/translations using RFdiffusion and PELE. 

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

This repository allows running the protein-protein PELE equilibration-production protocol and the DiffPELE protocol, which is a pipeline that integrates the next workflow, given a PDB protein-protein system:
- RFdiffusion (partial diffusion protocol) to generate an ensemble of poses with different interface backbone conformations.
- FastRelax to repack and minimize side chains in diffusion models.
- PELE to refine (adjust the binding mode and mainly side chains conformations) the repacked and minimized diffusion models.

### System preparation

Given your protein-protein system, keep only amino acids from chains of interest and remove waters, small molecules, etc.
Use the script `./DiffPELE/src/0_preprocess_system.py` to fastly do certain PDB operations (remove/save chains, remove hetatms, add TERs, renumber antibody insertion numbers...) although it is not compulsory. The prepared PDB must not have residue numbers gaps (1, 2, 3, 4, 6) or antibody insertion codes (1, 1A, 1B, 2...).

Finally, run this Schrodinger pipeline:
```bash
$SCHRODINGER/utilities/prepwizard ./DiffPELE/examples/4POU.pdb ./DiffPELE/examples/4POU_prep.pdb -rehtreat -disulfides -fillloops -fillsidechains -propka_pH 7.4 -minimize_adj_h -f OPLS_2005
```
where $SCHRODINGER stands for the path of the Maestro Schrodinger molecular modeling suite.

### Run DiffPELE in MNV

To run DiffPELE in MNV use our runner and adapt parameters to your input PDB and chains or input from terminal. See these 2 examples:
```bash
sbatch ./DiffPELE/diffpele_MNV.runner ./DiffPELE/examples/4POU.pdb B A
sbatch ./DiffPELE/diffpele_MNV.runner ./DiffPELE/examples/5C7X.pdb H,L A
```
Finally, check energy profiles in reports located at: `(system)_diffpele/(system)_diffusion_pele/outPROD/0`

### Run PELE in MNV

Similarly, to run PELE in MNV on protein-protein systems, use our runner and adapt parameters to your input PDB and chains or input from terminal. See these 2 examples:
```bash
sbatch ./DiffPELE/pele_MNV.runner ./DiffPELE/examples/4POU.pdb B A
sbatch ./DiffPELE/pele_MNV.runner ./DiffPELE/examples/5C7X.pdb H,L A
```
Finally, check energy profiles in reports located at: `(system)_peleBaseline/outPROD/0`
In this context, run PELE on a system with non-optimized interface conformations and on a system with optimized interface conformations to account for the energy baseline.

### Run DockQ
We use DockQ to measure (according to structure metrics) how unoptimized is a system in comparison to its optimized version. If there is one receptor chain and one ligand chain (thus one interface) DockQ can be run like that.

```bash
DockQ <a> <b>
```
where 'a' stands for the system with non-optimized interactions and 'b' stands for the optimized/native one. In case the receptor has 2 chains, you need to merge the receptor chains. See for further information: https://github.com/bjornwallner/DockQ/issues/33

To run DockQ on the best PELE poses (in binding energy), run our script:
```bash
bash ./DiffPELE/bin/evaluate_models.sh /path/to/PELE/dir /path/to/native_system.pdb receptorChain(s) ligandChain
bash ./DiffPELE/bin/evaluate_models.sh /path/to/PELE/dir /path/to/native_system.pdb B A
bash ./DiffPELE/bin/evaluate_models.sh /path/to/PELE/dir /path/to/native_system.pdb H,L A
```

### Run several DiffPELE iterations
Remember that in DiffPELE, as the last stage is running PELE: its trajectories or models have residue names that are different from the REF15 FF used in RFdiffusion. Thus, models after one DiffPELE epoch will need residue renaming before entering on a second epoch. Use our script:
```bash
python ./DiffPELE/src/5_postprocess_model.py /path/to/PDB/or/dirWithPDBs
```