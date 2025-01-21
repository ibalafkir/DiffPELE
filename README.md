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

My Script
Schrodinger

### Energy base line
Run PELE

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