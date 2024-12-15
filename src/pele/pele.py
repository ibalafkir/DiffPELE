import os 
import copy as cp
from src.utils.pdb_utils import PdbDf, PdbHandler
from src.utils.logger_factory import LoggerFactory
from src.utils.os_utils import display_full_dataframe
from src.utils.file_handling import write_yaml_data
import argparse
import textwrap
import shutil


logger = LoggerFactory.get_logger(__name__, "INFO")
display_full_dataframe()


def validate_and_get_chains(receptor_chains, ligand_chain):
    
    """
    Validate and get the chains to keep in the PDB file.
    
    Parameters:
    receptor_chains (str): Chain ID for the receptor protein chains.
    ligand_chain (str): Chain ID for the ligand protein chain.
    
    Returns:
    chains_to_keep (list): List of chains to keep in the PDB file
    """
    
    # Validate input chains
    
    if len(ligand_chain) > 1:
        raise ValueError("Only one ligand chain is allowed.")
    
    if len(receptor_chains) == 0:
        raise ValueError("At least one receptor chain must be inputted.")
    
    if len(ligand_chain) == 0:
        raise ValueError("At least one ligand chain must be inputted.")
    
    if len(receptor_chains) == 2 and ',' not in receptor_chains:
            raise ValueError("If two receptor chains are inputted, separate them with a comma.")
        
    if len(receptor_chains) > 3 and ',' in receptor_chains:
        raise ValueError("Only two receptor chains are allowed.")
    
    # Get chains to keep
    
    chains_to_keep = []
    
    if len(receptor_chains) == 1:
        chains_to_keep.append(receptor_chains)
        chains_to_keep.append(ligand_chain)
    else: # len(receptor_chains) == 3
        receptor_chains = receptor_chains.split(',')
        chains_to_keep.append(receptor_chains[0])
        chains_to_keep.append(receptor_chains[1])
        chains_to_keep.append(ligand_chain)
    
    print(chains_to_keep, receptor_chains, ligand_chain)
    return chains_to_keep, receptor_chains, ligand_chain


class PeleBuilder:
    
    def __init__(
        self,
        pdb_path: str,
        receptor_chains: str,
        ligand_chain: str,
        distance_cutOff: float = 12.0
    ):
        """
        Initialize the PeleBuilder class.
        
        Parameters:
        pdb_path (str): Path to the PDB file.
        receptor_chains (str): Chain ID for the receptor protein chains or chain.
        ligand_chain (str): Chain ID for the ligand protein chain.
        distance_cutOff (float): Distance cut-off for interaction analysis. Default is 12.0.
        """
        
        # Validate inputs
        self._validate_inputs(
            pdb_path=pdb_path,
            receptor_chains=receptor_chains,
            ligand_chain=ligand_chain,
            distance_cutOff=distance_cutOff
        )
        
        self.pdb_path = pdb_path
        self.receptor_chains = receptor_chains
        self.ligand_chain = ligand_chain
        self.distance_cutOff = distance_cutOff
        
        self.pdb_code = None
        self.pdb_file = None
        self.base_dir = None
        self.base_dir_pdb = None
        self.base_dir_control = None
        self.receptor_chain_1, self.receptor_chain_2 = None, None # in case there are 2 receptor chains
        
        self._get_pdb_info()
        self._organize_dir()
        
        
    def _validate_inputs(
        self,
        pdb_path: str,
        receptor_chains: str,
        ligand_chain: str,
        distance_cutOff: float
    ):
        
        # Validate string inputs
        for var_name, var_value in [
            ("pdb_path", pdb_path),
            ("receptor_chains", receptor_chains),
            ("ligand_chain", ligand_chain)
        ]:
            if not isinstance(var_value, str):
                raise TypeError(f"{var_name} should be a string")
        
        # Validate float inputs
        if not isinstance(distance_cutOff, float):
            raise TypeError(f"{distance_cutOff} should be a float")
        
         # Validate PDB
        if not os.path.exists(pdb_path):
            raise FileNotFoundError(f"File not found: {pdb_path}"
                                    )       
        # Validate distance cut-off
        if distance_cutOff <= 0:
            raise ValueError(f"{distance_cutOff} should be positive")
        
        # Validate chains in PDB
        atom_df = PdbDf(pdb_path)
        atom_df.get_atoms()
        atom_df_chains = atom_df.get_chains_id()
        if receptor_chains[0] not in atom_df_chains or receptor_chains[::-1][0] not in atom_df_chains: # since recepotr_chain can be A,B, verifiy if letters in atom_df
            raise ValueError(f"{receptor_chains} is not a valid chain ID in {pdb_path}")
        if ligand_chain not in atom_df_chains:
            raise ValueError(f"{ligand_chain} is not a valid chain ID in {pdb_path}")
        if receptor_chains == ligand_chain or ligand_chain in receptor_chains:
            raise ValueError(f"Chain IDs should be different")
        if receptor_chains == ' ' or ligand_chain == ' ':
            raise ValueError(f"Chain IDs should not be empty")
        if len(receptor_chains) == 2 or len(receptor_chains) > 3 or len(ligand_chain) > 1:
            raise ValueError(f"Chain IDs should be 1 character for the ligand chain and 1 (letter) \
                            or 3 (letter,letter) characters for the receptor chain")


    def _get_pdb_info(self):
        
        self.pdb_code = os.path.basename(self.pdb_path).split(".")[0]
        self.pdb_file = os.path.basename(self.pdb_path)

        if len(self.receptor_chains) > 1:
            self.receptor_chain_1 = self.receptor_chains[0]
            self.receptor_chain_2 = self.receptor_chains[::-1][0]
            logger.info(f"Detected two receptor chains: {self.receptor_chain_1} and {self.receptor_chain_2} and one ligand chain: {self.ligand_chain}")
        
        else:
            logger.info(f"Detected one receptor chain: {self.receptor_chains} and one ligand chain: {self.ligand_chain}")
            
        
    def _organize_dir(self):
        
        self.pdb_path = os.path.abspath(self.pdb_path)
        self.base_dir = self.pdb_path[:-4]+"_pele"
        self.base_dir_pdb = self.base_dir+"/pdbs"
        self.base_dir_control = self.base_dir+"/ctrl"
        print(self.base_dir_control)
        os.makedirs(self.base_dir, exist_ok=True)
        os.makedirs(self.base_dir_pdb, exist_ok=True)
        os.makedirs(self.base_dir_control, exist_ok=True)
        shutil.copy(self.pdb_path, self.base_dir_pdb)
        

    def create_nbdsuite_runner(self):

        content = f"""\
            #!/bin/bash
            #SBATCH -J {self.pdb_code}_sui
            #SBATCH --output=%x.out
            #SBATCH --error=%x.err
            #SBATCH --ntasks=1
            #SBATCH --qos=gp_debug
            #SBATCH --time=0:10:00
            #SBATCH --account=bsc72

            ml intel/2023.2.0
            ml cmake/3.25.1
            ml impi/2021.10.0
            ml mkl/2023.2.0
            ml miniconda/24.1.2
            ml boost/1.77.0-gcc

            source activate /gpfs/projects/bsc72/conda_envs/nbdsuite/0.5.0
            umask u=rwx,g=rwx,o=

            python -m nbdsuite.main input.yaml
        """
        
        content = textwrap.dedent(content)
        
        output_path = os.path.join(self.base_dir, "nbdsuite.sh")
        
        with open(output_path, "w") as file:
            file.write(content)
            
            
    def create_nbdsuite_input(self):
        
        content = f"""\
            system_data: pdbs/*.pdb
            cpus: 1
            name: processed
            pipeline:
                - block: topology_extractor
                - block: pele_pdb_preprocessor
                options:
                    set_unique_pdb_atom_names: True
                    fix_side_chains: True
                - block: topology_retriever
                options:
                    restore_original_names: True
        """

        content = textwrap.dedent(content)
        output_path = os.path.join(self.base_dir, "nbdsuite_input.yaml")

        with open(output_path, "w") as file:
            file.write(content)

    
    def create_equilibration_runner(self):
        
        content = f"""\
            #!/bin/bash
            #SBATCH --job-name {self.pdb_code[:2]}_eq
            #SBATCH --output=%x.out
            #SBATCH --error=%x.err
            #SBATCH --ntasks=16
            #SBATCH --time=2:00:00            
            #SBATCH -D .
            #SBATCH --cpus-per-task=1
            #SBATCH --qos=gp_debug
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

            python -m AdaptivePELE.adaptiveSampling control_files/adaptiveEQ.conf
        """
        
        content = textwrap.dedent(content)
        out_path = os.path.join(self.base_dir, "runEq.sh")
        
        with open(out_path, "w") as file:
            file.write(content)


    def create_production_runner(self):

        content = f"""\
                #!/bin/bash
                #SBATCH --job-name {self.pdb_code[:2]}_prod
                #SBATCH --output=%x.out
                #SBATCH --error=%x.err
                #SBATCH --ntasks=64
                #SBATCH --time=24:00:00            
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

                python -m AdaptivePELE.adaptiveSampling control_files/adaptivePROD.conf
        """

        content = textwrap.dedent(content)
        output_path = os.path.join(self.base_dir, "runProd.sh")

        with open(output_path, "w") as file:
            file.write(content)

    
    def create_control_adaptive_equilibration(self): # #"""+self.pdb_code+"""_eq",
                # TODO
        content = """\
            {{
                "generalParams" : {{
                    "restart": true,
                    "outputPath": "{self.pdb_code}_eq",
                    "initialStructures" : ["processed/3_topology_retriever/systems/*_1.pdb"]
                }},
                "spawning" : {{
                    "type" : "independent",
                    "params" : {{
                        "reportFilename" : "report",
                        "metricColumnInReport" : 5,
                        "epsilon": 0.50,
                        "T":1000
                    }},
                    "density" :{{
                        "type": "null"
                    }}
                }},
                "simulation": {{
                    "type" : "pele",
                    "params" : {{
                        "iterations" : 1,
                        "peleSteps" : 40,
                        "processors" : 16,
                        "runEquilibration" : false,
                        "equilibrationLength" : 1,
                        "seed": 12345,
                        "executable": "/gpfs/projects/bsc72/PELE++/mnv/1.8.0/bin/PELE-1.8_mpi",
                        "data": "/gpfs/projects/bsc72/PELE++/mnv/1.8.0/Data",
                        "documents": "/gpfs/projects/bsc72/PELE++/mnv/1.8.0/Documents",
                        "useSrun": true,
                        "controlFile" : "control_files/peleEQ.conf"
                    }}
                }},
                "clustering" : {{
                    "type" : "rmsd",
                    "params" : {{
                        "ligandChain" : "?",
                        "alternativeStructure" : true,
                        "contactThresholdDistance" : 8
                    }},
                    "thresholdCalculator":{{
                        "type" : "heaviside",
                        "params" : {{
                        "values" : [2.0, 3, 5],
                        "conditions": [0.20, 0.12, 0]
                        }}
                    }}
                }}
            }}
        """
        
        content = textwrap.dedent(content)
        output_path = os.path.join(self.base_dir_control, "adEq.conf")

        with open(output_path, 'w') as file:
            file.write(content)

    def create_control_adaptive_production(self):
        
        content = """\
            {
                "generalParams" : {
                    "restart": true,
                    "outputPath":" """+self.pdb_code+"""_eq",
                "initialStructures" : ["pdbs/eq.pdb"]    },
                "spawning" : {
                    "type" : "independent",
                    "params" : {
                        "reportFilename" : "report",
                        "metricColumnInReport" : 5,
                        "epsilon": 0.50,
                        "T":1000
                },
                    "density" :{
                        "type": "null"
                    }
                },
                "simulation": {
                    "type" : "pele",
                    "params" : {
                        "iterations" : 1,
                        "peleSteps" : 200,
                        "processors" : 64,
                        "runEquilibration" : false,
                        "equilibrationLength" : 1,
                        "seed": 12345,
                        "executable": "/gpfs/projects/bsc72/PELE++/mnv/1.8.0/bin/PELE-1.8_mpi",
                        "data": "/gpfs/projects/bsc72/PELE++/mnv/1.8.0/Data",
                        "documents": "/gpfs/projects/bsc72/PELE++/mnv/1.8.0/Documents",
                        "useSrun": true,
                        "controlFile" : "control_files/pelePROD.conf"
                    }
                },
                "clustering" : {
                    "type" : "rmsd",
                    "params" : {
                        "ligandChain" : " """+self.ligand_chain+""" ",
                        "alternativeStructure" : true,
                        "contactThresholdDistance" : 8
                    },
                    "thresholdCalculator":{
                        "type" : "heaviside",
                        "params" : {
                        "values" : [2.0, 3, 5],
                        "conditions": [0.20, 0.12, 0]
                        }
                    }
                }
            }     
        """
        
        content = textwrap.dedent(content)
        output_path = os.path.join(self.base_dir_control, "adProd.conf")
        
        with open(output_path, 'w') as file:
            file.write(content)

if __name__ == "__main__":
    
    test = PeleBuilder(
        pdb_path="/home/ibalakfi/Desktop/4pou.pdb",
        receptor_chains="A",
        ligand_chain="B",
        distance_cutOff=12.0
    )
    
    
    
    
    
    test.create_nbdsuite_input()
    test.create_nbdsuite_runner()
    test.create_equilibration_runner()
    test.create_production_runner()
    test.create_control_adaptive_equilibration()
    
    # TODO comment the whole script