import os 
import copy as cp
from src.utils.pdb_utils import PdbDf, PdbHandler
from src.utils.logger_factory import LoggerFactory
from src.utils.os_utils import display_full_dataframe
from src.utils.file_handling import write_yaml_data
from src.interface import InterfaceAnalyzer
import argparse
import textwrap
import shutil
import json
import glob


logger = LoggerFactory.get_logger(__name__, "INFO")
display_full_dataframe()

#####
### TODO:
# Completar el create pele conf prod
# en adaptive prod, si el modo es diffusion, que coja los pdb equilibrados de pdbsEQ que creara el script de equilibrar pele
####3

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
    
    return chains_to_keep, receptor_chains, ligand_chain


class PeleBuilder:
    
    def __init__(
        self,
        pdb_path: str,
        receptor_chains: str,
        ligand_chain: str,
        distance_cutOff: float = 12.0,
        mode: str = "single",
        diffmodels_path: str = None
    ):
        """
        Initialize the PeleBuilder class.
        Build PELE configuration and control files for a given PDB file,
        the receptor and ligand chains, and a distance cut-off.
        Resulting files will be run in MareNostrumV, the PELE++ and 
        AdaptivePELE conda environments
        
        Parameters:
        pdb_path (str): Path to the PDB file.
        receptor_chains (str): Chain ID for the receptor protein chains or chain.
        ligand_chain (str): Chain ID for the ligand protein chain.
        distance_cutOff (float): Distance cut-off for interaction analysis. Default is 12.0.
        mode (str): Mode to run PELE. 
                    "single" to run for a single system. Default
                    "diffusion" to run for multiple diffusion models originated from a single system.
                    in this case, interface metrics are obtained from this first system and thus its
                    path needs to be provided.
        diffmodels_path (str): Path to the diffusion models. Required only if mode is "diffusion".
                               It must contain 40 PDB diffusion models.
        """
        
        # Validate inputs
        self._validate_inputs(
            pdb_path=pdb_path,
            receptor_chains=receptor_chains,
            ligand_chain=ligand_chain,
            distance_cutOff=distance_cutOff,
            mode=mode,
        )
        
        # Set attributes
        self.pdb_path = pdb_path
        self.receptor_chains = receptor_chains
        self.ligand_chain = ligand_chain
        self.distance_cutOff = distance_cutOff
        self.mode = mode
        self.diffmodels_path = diffmodels_path
        
        # Set attributes for use in other methods
        self.pdb_code = None
        self.pdb_file = None
        self.base_dir = None
        self.base_dir_pdb = None
        self.base_dir_control = None
        self.receptor_chain_1, self.receptor_chain_2 = None, None # in case there are 2 receptor chains
        
        # Get PDB info and make directory structure
        self._get_pdb_info()
        self._organize_dir()
        
        
    def _validate_inputs(
        self,
        pdb_path: str,
        receptor_chains: str,
        ligand_chain: str,
        distance_cutOff: float,
        mode: str,
        diffmodels_path: str = None
    ):
        
        # Validate string inputs
        for var_name, var_value in [
            ("pdb_path", pdb_path),
            ("receptor_chains", receptor_chains),
            ("ligand_chain", ligand_chain),
            ("mode", mode)
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
        
        # Validate mode
        if mode not in ["single", "diffusion"]:
            raise ValueError(f"{mode} is not a valid mode. Choose 'single' or 'diffusion'")
        
        if mode == "single" and diffmodels_path is not None:
            raise ValueError("diff_path_models is not required for mode 'single' since you \
                              are running a single system")       
        if mode == "diffusion" and diffmodels_path is None:
            raise ValueError("diff_path_models is required for mode 'diffusion'")
                
        # Validate diffmodels_path
        if diffmodels_path != None and not isinstance(diffmodels_path, str):
            raise TypeError(f"{diffmodels_path} must be a string")
        if mode == "diffusion":
            if not os.path.exists(diffmodels_path):
                raise FileNotFoundError(f"File not found: {diffmodels_path}, required for mode 'diffusion'")
            diffmodels_path_files = glob.glob(os.path.join(diffmodels_path, "*.pdb"))
            if len(diffmodels_path_files) == 0:
                raise FileNotFoundError(f"No PDB files found in {diffmodels_path}")
            
            if len(diffmodels_path_files) > 0 and len(diffmodels_path_files) != 2:
                raise ValueError(f"The current diffusion mode expects to run with 40 systems.") # TODO change in future
        
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
        
        # Get PDB 4 letter code and file name
        self.pdb_code = os.path.basename(self.pdb_path).split(".")[0]
        self.pdb_file = os.path.basename(self.pdb_path)

        # Get chains info
        if len(self.receptor_chains) > 1:
            self.receptor_chain_1 = self.receptor_chains[0]
            self.receptor_chain_2 = self.receptor_chains[::-1][0]
            logger.info(f"Detected two receptor chains: {self.receptor_chain_1} and {self.receptor_chain_2}; and one ligand chain: {self.ligand_chain}")
            self.receptor_chains = [self.receptor_chain_1, self.receptor_chain_2]
            logger.info(f"Receptor chains: {self.receptor_chains}")
        else:
            logger.info(f"Detected one receptor chain: {self.receptor_chains}; and one ligand chain: {self.ligand_chain}")
            
        
    def _organize_dir(self):
    
        ### Create base directory and subdirectories ###
        
        # Absolute path to PDB
        self.pdb_path = os.path.abspath(self.pdb_path)
        
        # Base PELE directory
        self.base_dir = self.pdb_path[:-4]+"_pele"
        
        # Subdirectories in base PELE directory
        self.base_dir_pdb = self.base_dir+"/pdbs"
        self.base_dir_control = self.base_dir+"/ctrl"
        
        # Make directories
        os.makedirs(self.base_dir, exist_ok=True)
        os.makedirs(self.base_dir_pdb, exist_ok=True)
        os.makedirs(self.base_dir_control, exist_ok=True)
        
        # Copy input PDB file to base_dir/pdbs (PELE runners will look for PDBs there)
        if self.mode == "single":
            shutil.copy(self.pdb_path, self.base_dir_pdb)
        # If mode is diffusion, copy all PDBs in diffmodel_path to base_dir/pdbs
        else: 
            pdb_files = glob.glob(os.path.join(self.diffmodels_path, "*.pdb"))
            for pdb_file in pdb_files:
                shutil.copy(pdb_file, self.base_dir_pdb)
                

    def create_nbdsuite_runner(self):

        # Content for runner
        content = f"""\
            #!/bin/bash
            #SBATCH -J {self.pdb_code}_sui
            #SBATCH --output=%x.out
            #SBATCH --error=%x.err
            #SBATCH --ntasks=____noOfCPUs____
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

        # Adjustments  
        if self.mode == "single":
            content = content.replace("____eqJobName____", f"{self.pdb_code[:2]}_sui")
            content = content.replace("____noOfCPUs____", "1")    
        else:
            content = content.replace("____eqJobName____", f"{self.pdb_code[:2]}_df_sui")
            content = content.replace("____noOfCPUs____", "40")
        
        # Remove leading spaces and write content into file
        content = textwrap.dedent(content)
        output_path = os.path.join(self.base_dir, "nbdsuite.sh")
        with open(output_path, "w") as file:
            file.write(content)
            
            
    def create_nbdsuite_input(self):
        
        # Content for input.yaml
        content = f"""\
            system_data: pdbs/*.pdb
            cpus: ____noOfCPUs____
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

        # Adjustments  
        if self.mode == "single":
            content = content.replace("____noOfCPUs____", "1")
        else: # mode == "diffusion"
            content = content.replace("____noOfCPUs____", "40")

        # Remove leading spaces and write content into file
        content = textwrap.dedent(content)
        output_path = os.path.join(self.base_dir, "nbdsuite_input.yaml")
        with open(output_path, "w") as file:
            file.write(content)

    
    def create_equilibration_runner(self):
        
        # Content for runner
        content = f"""\
            #!/bin/bash
            #SBATCH --job-name ____eqJobName____
            #SBATCH --output=%x.out
            #SBATCH --error=%x.err
            #SBATCH --ntasks=____noOfCPUs____
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

            python -m AdaptivePELE.adaptiveSampling ctrl/adEQ.conf
        """
        
        # Adjustments
        if self.mode == "single":
            content = content.replace("____eqJobName____", f"{self.pdb_code[:2]}_eq")
            content = content.replace("____noOfCPUs____", "16")

        else: # mode == "diffusion"
            content = content.replace("____eqJobName____", f"{self.pdb_code[:2]}_df_eq")
            content = content.replace("____noOfCPUs____", "201")
            
        # Remove leading spaces and write content into file
        content = textwrap.dedent(content)
        out_path = os.path.join(self.base_dir, "runEq.sh")
        with open(out_path, "w") as file:
            file.write(content)


    def create_production_runner(self):

        # Content for runner
        content = f"""\
                #!/bin/bash
                #SBATCH --job-name ____prodJobName____
                #SBATCH --output=%x.out
                #SBATCH --error=%x.err
                #SBATCH --ntasks=____noOfCPUs____
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

                python -m AdaptivePELE.adaptiveSampling ctrl/adProd.conf
        """
        
        # Adjustments  
        if self.mode == "single":
            content = content.replace("____prodJobName____", f"{self.pdb_code[:2]}_prod")
            content = content.replace("____noOfCPUs____", "64")
        else: # mode == "diffusion"
            content = content.replace("____prodJobName____", f"{self.pdb_code[:2]}_df_prod")
            content = content.replace("____noOfCPUs____", "209")
        
        # Remove leading spaces and write content into file
        content = textwrap.dedent(content)
        output_path = os.path.join(self.base_dir, "runProd.sh")
        with open(output_path, "w") as file:
            file.write(content)


    def create_control_adaptive_equilibration(self):

        # Content for control file
        content = """\
        {
            "generalParams" : {
                "restart": true,
                "outputPath": "____outputPathNameToChange____",
            "initialStructures" : ["processed/3_topology_retriever/systems/*_1.pdb"]    },
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
                    "iterations" : ____noOfIter____,
                    "peleSteps" : ____noOfPeleSteps____,
                    "processors" : ____noOfCPUs____,
                    "runEquilibration" : false,
                    "equilibrationLength" : 1,
                    "seed": 12345,
                    "executable": "/gpfs/projects/bsc72/PELE++/mnv/1.8.0/bin/PELE-1.8_mpi",
                    "data": "/gpfs/projects/bsc72/PELE++/mnv/1.8.0/Data",
                    "documents": "/gpfs/projects/bsc72/PELE++/mnv/1.8.0/Documents",
                    "useSrun": true,
                    "controlFile" : "control_files/peleEQ.conf"
                }
            },
            "clustering" : {
                "type" : "rmsd",
                "params" : {
                    "ligandChain" : "____ligandChainToChange____",
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
        
        # Adjustments   
        if self.mode == "single":
            content = content.replace("____outputPathNameToChange____", f"{self.pdb_code}_eq")
            content = content.replace("____noOfIter____", "1")
            content = content.replace("____noOfPeleSteps____", "40")
            content = content.replace("____noOfCPUs____", "16")            
            content = content.replace("____ligandChainToChange____", self.ligand_chain)
        else: # mode == "diffusion"
            content = content.replace("____outputPathNameToChange____", f"{self.pdb_code}_df_eq")
            content = content.replace("____noOfIter____", "1")
            content = content.replace("____noOfPeleSteps____", "40")
            content = content.replace("____noOfCPUs____", "201")

        # Remove leading spaces and write content into file
        content = textwrap.dedent(content)
        output_path = os.path.join(self.base_dir_control, "adEq.conf")
        with open(output_path, 'w') as file:
            file.write(content)


    def create_control_adaptive_production(self):

        # Content for control file
        content = """\
            {
                "generalParams" : {
                    "restart": true,
                    "outputPath":"____outputPathNameToChange____",
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
                        "iterations" : ____noOfIter____,
                        "peleSteps" : ____noOfPeleSteps____,
                        "processors" : ____noOfCPUs____,
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
                        "ligandChain" : "____ligandChainToChange____",
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
        
        # Adjustments

        if self.mode == "single":
            content = content.replace("____outputPathNameToChange____", f"{self.pdb_code}_eq")            
            content = content.replace("____noOfIter____", "1")
            content = content.replace("____noOfPeleSteps____", "200")
            content = content.replace("____noOfCPUs____", "64")
            content = content.replace("____ligandChainToChange____", self.ligand_chain)
        else: # mode == "diffusion"
            content = content.replace("____outputPathNameToChange____", f"{self.pdb_code}_df_eq")
            content = content.replace("____noOfIter____", "1")
            content = content.replace("____noOfPeleSteps____", "200")
            content = content.replace("____noOfCPUs____", "209")
            content = content.replace("____ligandChainToChange____", self.ligand_chain)

        # Remove leading spaces and write content into file
        content = textwrap.dedent(content)
        output_path = os.path.join(self.base_dir_control, "adProd.conf")
        with open(output_path, 'w') as file:
            file.write(content)


    def create_pele_conf_equilibration(self):
        
        # Content for control file
        content = """\
        {
        "licenseDirectoryPath" : "/gpfs/projects/bsc72/PELE++/license",
        "simulationLogPath" : "$OUTPUT_PATH/logFile.txt",  
        "Initialization" : {
            "allowMissingTerminals": true,
            "ForceField" : "OPLS2005",
            "MultipleComplex": [ $COMPLEXES ],
            "Solvent" : {
                "ionicStrength" : 0.15, "solventType" : "VDGBNP", "useDebyeLength" : true }
        },
        "verboseMode": true,
        "commands" : [
            {
                "commandType" : "peleSimulation",
                "RandomGenerator" : { "seed" : $SEED },
                "selectionToPerturb" : { "chains" : { "names" : [ "____ligandChainToChange____" ] } },
                "PELE_Output" : {
                    "savingFrequencyForAcceptedSteps" : 1,
                    "savingMode" : "savingTrajectory",
                    "reportPath": "$OUTPUT_PATH/report",
                    "trajectoryPath": "$OUTPUT_PATH/trajectory.pdb"
                },
                "PELE_Parameters" : {
                    "anmFrequency" : 4,
                    "sideChainPredictionFrequency" : 1,
                    "minimizationFrequency" : 1,
                    "activateProximityDetection": false,
                    "temperature": 1500,
                    "numberOfPeleSteps": $PELE_STEPS
                },
                "Perturbation": {
                        "perturbationType":"naive",
                        "translationDirection": "random",
                        "rotationAngles": "nonCoupled",
                        "parameters": {
                            "peleRegionType": "interfaceLinks",
                            "steeringUpdateFrequency": 1,
                            "influenceRange": 3, 
                    "perturbAllAtOnce": true
                        }   
                    },
                "ANM" : {
                    "algorithm": "CARTESIANS", "nodes": { "atoms": { "names": [ "_CA_" ]} },
                    "ANMMinimizer" : {
                    "algorithm" : "TruncatedNewton",
                    "parameters" : {
                        "MaximumMinimizationIterations" : 1,
                        "MaximumNewtonIterations" : 20,
                        "MinimumRMS" : 1.0,
                        "alphaUpdated" : false,
                        "nonBondingListUpdatedEachMinStep" : false
                    }
                    },
                    "options" : {
                    "directionGeneration" : "random",
                    "modesMixingOption" : "mixMainModeWithOthersModes",
                    "pickingCase" : "RANDOM_MODE"
                    },
                    "parameters" : {
                    "displacementFactor" : 0.25,
                    "eigenUpdateFrequency" : 1000000,
                    "mainModeWeightForMixModes" : 0.5,
                    "modesChangeFrequency" : 4,
                    "numberOfModes": 6,
                    "relaxationSpringConstant" : 0.0
                    }
            },

                "SideChainPrediction" : {
                    "algorithm" : "zhexin",
                    "parameters" : { "discardHighEnergySolutions" : false, "resolution": 10, "randomize" : false, "numberOfIterations": 1 }
                },

                "Minimizer" : {
                    "algorithm" : "TruncatedNewton",
                    "parameters" : { "MaximumMinimizationIterations" : 20, "MaximumNewtonIterations" : 1, "MinimumRMS" : 0.1, "alphaUpdated" : false, "nonBondingListUpdatedEachMinStep" : false }
                },
                "PeleTasks" : [
                    {

                        "metrics" : [

                                { "type": "bindingEnergy",
                                    "tag": "Binding_energy",
                                    "boundPartSelection": { "chains": { "names": ["____ligandChainToChange____"] } },
                                    "allowMultipleBindingSelection" : true

                                },

                                { "type": "sasa",
                                    "tag": "SASA_ligand",
                                    "selection": { "chains": { "names": ["____ligandChainToChange____"] } }

                                },

                                {
                                "type":"com_distance",
                                "tag":"ParaEpi_distance",
                                "selection_group_1":{
                                    "links": { "ids":____receptorInteractingRes____}},
                                },
                                "selection_group_2":{
                                    "links": { "ids":____ligandInteractingRes____}},
                                },

                                {
                                "type": "rmsd",
                                "tag": "L-RMSD",
                                "Native": {
                                        "path": "./pdbs/____pdbNameToChange____",},
                                "selection": {
                                        "chains": {"names": [ "____ligandChainToChange____" ] },
                                        "atoms": {"names": [ "_CA_" ]}},
                                "doSuperposition": true,
                                "superpositionSelection": {
                                        "chains": {"names": ____receptorChainsToChange____ }}
                                },
                                
                                { "tag" : "rand", "type" : "random" },
                                { "tag" : "rand2", "type" : "random" },
                                { "tag" : "rand1", "type" : "random" }
                            ] ,


                    "parametersChanges" : [
                        { "ifAnyIsTrue": [ "rand1 >= 0.5" ],
                        "doThesechanges": { "Perturbation::parameters": {  "rotationScalingFactor": 0.01 } },
                        "otherwise": { "Perturbation::parameters": { "rotationScalingFactor": 0.01 } }
                        },
                        { "ifAnyIsTrue": [ "rand >= 0.5" ],
                        "doThesechanges": { "Perturbation::parameters": { "translationRange": 0.02, "numberOfTrials" : 2, "numberOfStericTrials": 500  } },
                        "otherwise": { "Perturbation::parameters": { "translationRange": 0.02, "numberOfTrials" : 2, "numberOfStericTrials": 500 } }
                        }
                    ]
                    }
                ]
                }
                ]
        }
        """

        # Adjustments
        
        if self.mode == "single":
        
            # Get interacting residues
            # only interface between receptor chain 1 and ligand chain in case the receptor has 2 chains
            if len(self.receptor_chains) == 1 and self.receptor_chain_1 == None and self.receptor_chain_2 == None:
                analyzer = InterfaceAnalyzer(
                    pdb_path=self.pdb_path,
                    receptor_chain=self.receptor_chains,
                    ligand_chain=self.ligand_chain,
                )
            else:
                analyzer = InterfaceAnalyzer(
                    pdb_path=self.pdb_path,
                    receptor_chain=self.receptor_chain_1,
                    ligand_chain=self.ligand_chain,
                )
            interaction_matrix = analyzer.get_interaction_matrix()
            interaction_matrix_receptor_chain_expanded, interaction_matrix_ligand_chain_expanded \
                = analyzer.expand_interface(neighborhood=3)
            interaction_matrix_receptor_chain_expanded_lst, interaction_matrix_ligand_chain_expanded_lst \
                = analyzer.get_interacting_residues(mode='e')
            interaction_matrix_receptor_chain_expanded_lst_pele = \
                [f"{item[0]}:{item[1:]}" for item 
                in interaction_matrix_receptor_chain_expanded_lst]
            interaction_matrix_ligand_chain_expanded_lst_pele = \
                [f"{item[0]}:{item[1:]}" for item 
                in interaction_matrix_ligand_chain_expanded_lst]
                
            # Adjustments
            # Json.dumps serves to represent python lists ['A', 'B'...] as json lists ["A", "B"...]
            content = content.replace("____ligandChainToChange____", self.ligand_chain)
            content = content.replace("____receptorInteractingRes____", json.dumps(interaction_matrix_receptor_chain_expanded_lst_pele))
            content = content.replace("____ligandInteractingRes____", json.dumps(interaction_matrix_ligand_chain_expanded_lst_pele))
            content = content.replace("____pdbNameToChange____", self.pdb_file)
            content = content.replace("____receptorChainsToChange____", json.dumps(self.receptor_chains))
        
            # Remove leading spaces and write content into file
            content = textwrap.dedent(content)
            output_path = os.path.join(self.base_dir_control, "peEq.conf")
            with open(output_path, 'w') as file:
                file.write(content)

        else: # mode == "diffusion"
            
            print("log")


    def create_pele_conf_production(self):
        
        pass




if __name__ == "__main__":
    
    test = PeleBuilder(
        pdb_path="/home/ibalakfi/Desktop/5C7X.pdb",
        receptor_chains="H,L",
        ligand_chain="A",
        distance_cutOff=12.0,
        mode = "single"
    )
    
    
    
    
    
    #test.create_nbdsuite_input()
    #test.create_nbdsuite_runner()
    #test.create_equilibration_runner()
    #test.create_production_runner()
    #test.create_control_adaptive_equilibration()
    #test.create_control_adaptive_production()
    test.create_pele_conf_equilibration()
    
    # tests
    # TODO single mode for 4pou
    # TODO diffusion model for 4pou
    # TODO signel model for 5c7x
    # TODO diffusion model for 5c7x