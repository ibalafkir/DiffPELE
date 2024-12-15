# DiffPELE
This pipeline optimizes flexible protein-protein and antibody-antigen interactions. This is done through interface backbone diffusion and Monte Carlo rotations/translations using RFdiffusion and PELE. The objective is to accurately model the conformations of the regions involved in the antibody-antigen or protein-protein interactions.

## Installation (WIP)
### Environment manual setup
```bash
$ git clone https://github.com/ibalafkir/diffpele.git
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
