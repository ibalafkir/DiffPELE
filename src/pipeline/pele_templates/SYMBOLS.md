# PELE Template Symbols

The following symbols are used in the default configuration templates. These placeholders are replaced with user-specified values during execution.

| Symbol               | Description                                           | Context           |
|----------------------|-------------------------------------------------------|-------------------|
| `_SJRN_`             | **Nbdsuite Job Runner Name**: Specifies the job runner name for NBDSuite. | General           |
| `_EJRN_`             | **Equilibration Job Runner Name**: Specifies the job runner name for the equilibration part of PELE. | General           |
| `_nCPUs_`            | **Number of CPUs**: Indicates the number of CPUs to be used in equilibration or production simulations. | General           |
| `_PJRN_`             | **Production Job Runner Name**: Specifies the job runner name for the production stage. | General           |
| `_OPN_`              | **Output Path Name**: Specifies the output path in adaptivePELE. | Adaptive          |
| `_nE_`               | **Number of Epochs**: Indicates the total number of epochs/iterations. | General           |
| `_nS_`               | **Number of Steps**: Indicates the total number of steps. | Adaptive           |
| `_ligandChainName_`  | **Ligand Chain Name**: Specifies the name of the ligand chain. | Adaptive, Control |
| `_receptorChainsNames_` | **Receptor Chains Names**: Specifies the names of receptor chains. | Adaptive, Control |
| `_RIR_`              | **Receptor Interacting Residues**: Indicates the residues on the receptor involved in interactions. | General           |
| `_LIR_`              | **Ligand Interacting Residues**: Indicates the residues on the ligand involved in interactions. | General           |
| `_PDB_`              | **PDB File for RMSD**: Specifies the Protein Data Bank (PDB) file used for Root Mean Square Deviation (RMSD) calculations. | General           |

## Usage

The pipeline takes template files containing the symbols above and modifies them based on user-provided inputs.