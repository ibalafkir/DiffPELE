import os 
from src.utils.pdb_utils import PdbDf
from src.utils.logger_factory import LoggerFactory
from src.utils.os_utils import display_full_dataframe
from src.pipeline.interface import InterfaceAnalyzer
import textwrap
import shutil
import json


logger = LoggerFactory.get_logger(__name__, "INFO")
display_full_dataframe()


class PeleSetup:
    
    def __init__(
        self,
        pdb_path: str,
        receptor_chains: str,
        ligand_chain: str,
        distance_cutOff: float = 12.0
    ):
        """
        Initialize the PeleSetup class.
        Build PELE runners, configuration and control files for a given PDB file,
        the receptor and ligand chains, and a distance cut-off.
        Resulting files will be run in MareNostrumV, the PELE++ and 
        AdaptivePELE conda environments
        
        Parameters:
        pdb_path (str): Path to the PDB file. Necessary to compute the com_distance metric of PELE.
        receptor_chains (str): Chain ID for the receptor protein chains or chain.
        ligand_chain (str): Chain ID for the ligand protein chain.
        distance_cutOff (float): Distance cut-off for interaction analysis. Default is 12.0.
                                 Useful to compute the distance between the interacting residues.
        """
        
        # Validate inputs
        self._validate_inputs(
            pdb_path=pdb_path,
            receptor_chains=receptor_chains,
            ligand_chain=ligand_chain,
            distance_cutOff=distance_cutOff
        )
        
        # Set attributes
        self.pdb_path = pdb_path
        self.receptor_chains = receptor_chains
        self.ligand_chain = ligand_chain
        self.distance_cutOff = distance_cutOff
        
        # Define attributes for use in other methods
        self.pdb_code = None
        self.pdb_file = None
        self.base_dir = None
        self.base_dir_pdbs = None
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
        if receptor_chains[0] not in atom_df_chains or receptor_chains[::-1][0] not in atom_df_chains:
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
            # logger.info(f"Receptor chains: {self.receptor_chains}")
        else:
            logger.info(f"Detected one receptor chain: {self.receptor_chains}; and one ligand chain: {self.ligand_chain}")
            
        
    def _organize_dir(self):
    
        ### Create base directory and subdirectories ###
        
        # Absolute path to PDB
        self.pdb_path = os.path.abspath(self.pdb_path)
        
        # Base PELE directory
        self.base_dir = self.pdb_path[:-4]+"_pele"
        
        # Subdirectories in base PELE directory
        self.base_dir_pdbs = self.base_dir+"/pdbs"
        self.base_dir_control = self.base_dir+"/ctrl"
        self.base_dir_ref = self.base_dir+"/ref"
        self.base_dir_pdbsEQ = self.base_dir+"/pdbsEQ"
        
        # Make directories
        os.makedirs(self.base_dir) # terminate execution if already existing
        os.makedirs(self.base_dir_pdbs, exist_ok=True)
        os.makedirs(self.base_dir_control, exist_ok=True)
        os.makedirs(self.base_dir_ref, exist_ok=True)
        os.makedirs(self.base_dir_pdbsEQ, exist_ok=True)
        
        # Copy input PDB file to base_dir/pdbs (PELE runners will look for PDBs there)
        shutil.copy(self.pdb_path, self.base_dir_pdbs)
        shutil.copy(self.pdb_path, self.base_dir_ref)
                

    def create_nbdsuite_runner(self, suiteJobRunnerName, nCPUs):

        # Content for runner
        template_path = os.path.join(os.path.dirname(__file__), "pele_templates/nbdsuite.sh")
        with open(template_path, "r") as file:
            content = file.read()

        # Adjustments  
        content = content.replace("_SJRN_", suiteJobRunnerName)
        content = content.replace("_nCPUs_", str(nCPUs)) 
        
        # Remove leading spaces and write content into file
        content = textwrap.dedent(content)
        output_path = os.path.join(self.base_dir, "nbdsuite.sh")
        with open(output_path, "w") as file:
            file.write(content)
            
            
    def create_nbdsuite_input(self, nCPUs):
        
        # Content for input.yaml
        template_path = os.path.join(os.path.dirname(__file__), "pele_templates/nbdsuite_input.yaml")
        with open(template_path, "r") as file:
            content = file.read()

        # Adjustments  
        content = content.replace("_nCPUs_", str(nCPUs))

        # Remove leading spaces and write content into file
        content = textwrap.dedent(content)
        output_path = os.path.join(self.base_dir, "nbdsuite_input.yaml")
        with open(output_path, "w") as file:
            file.write(content)


    def create_equilibration_runner(self, nCPUs, equiJobRunnerName):
        
        # Content for runner
        content = os.path.join(os.path.dirname(__file__), "pele_templates/runEq.sh")
        with open(content, "r") as file:
            content = file.read()
        
        # Adjustments
        content = content.replace("_EJRN_", equiJobRunnerName)
        content = content.replace("_nCPUs_", str(nCPUs))
            
        # Remove leading spaces and write content into file
        content = textwrap.dedent(content)
        out_path = os.path.join(self.base_dir, "runEq.sh")
        with open(out_path, "w") as file:
            file.write(content)


    def create_production_runner(self, nCPUs, prodJobRunnerName):

        # Content for runner
        content = os.path.join(os.path.dirname(__file__), "pele_templates/runProd.sh")
        with open(content, "r") as file:
            content = file.read()
        
        # Adjustments  
        content = content.replace("_PJRN_", prodJobRunnerName)
        content = content.replace("_nCPUs_", str(nCPUs))

        # Remove leading spaces and write content into file
        content = textwrap.dedent(content)
        output_path = os.path.join(self.base_dir, "runProd.sh")
        with open(output_path, "w") as file:
            file.write(content)


    def create_control_adaptive_equilibration(self, outputPathName, nEpochs, nSteps, nCPUs):

        # Content for control file
        content = os.path.join(os.path.dirname(__file__), "pele_templates/adEq.conf")
        with open(content, "r") as file:
            content = file.read()
        
        # Adjustments   
        content = content.replace("_OPN_", outputPathName)
        content = content.replace("_nE_", str(nEpochs))
        content = content.replace("_nS_", str(nSteps))
        content = content.replace("_nCPUs_", str(nCPUs))            
        content = content.replace("_ligandChainName_", self.ligand_chain)

        # Remove leading spaces and write content into file
        content = textwrap.dedent(content)
        output_path = os.path.join(self.base_dir_control, "adEq.conf")
        with open(output_path, 'w') as file:
            file.write(content)


    def create_control_adaptive_production(self, outputPathName, nEpochs, nSteps, nCPUs):

        # Content for control file
        content = os.path.join(os.path.dirname(__file__), "pele_templates/adProd.conf")
        with open(content, "r") as file:
            content = file.read()
        
        # Adjustments   
        content = content.replace("_OPN_", outputPathName)
        content = content.replace("_nE_", str(nEpochs))
        content = content.replace("_nS_", str(nSteps))
        content = content.replace("_nCPUs_", str(nCPUs))            
        content = content.replace("_ligandChainName_", self.ligand_chain)

        # Remove leading spaces and write content into file
        content = textwrap.dedent(content)
        output_path = os.path.join(self.base_dir_control, "adProd.conf")
        with open(output_path, 'w') as file:
            file.write(content)


    def create_pele_conf_equilibration(self):
        
        # Content for control file
        content = os.path.join(os.path.dirname(__file__), "pele_templates/peEq.conf")
        with open(content, "r") as file:
            content = file.read()        

        # Get interacting residues
        # only interface between receptor chain 1 and ligand chain in case the receptor has 2 chains.
        # if the receptor is an antibody, input the heavy chain
        if len(self.receptor_chains) == 1 and self.receptor_chain_1 == None and self.receptor_chain_2 == None:
                analyzer = InterfaceAnalyzer(
                    pdb_path=self.pdb_path,
                    receptor_chain=self.receptor_chains,
                    ligand_chain=self.ligand_chain,
                )
        else: # receptor has 2 chains
            analyzer = InterfaceAnalyzer(
                pdb_path=self.pdb_path,
                receptor_chain=self.receptor_chain_1,
                ligand_chain=self.ligand_chain,
                )
                
        _ = analyzer.get_interaction_matrix()
        _ , _ = analyzer.expand_interface(neighborhood=3)
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
        content = content.replace("_ligandChainName_", self.ligand_chain)
        content = content.replace("_RIR_", json.dumps(interaction_matrix_receptor_chain_expanded_lst_pele))
        content = content.replace("_LIR_", json.dumps(interaction_matrix_ligand_chain_expanded_lst_pele))
        content = content.replace("_PDB_", self.pdb_file)
        
        if len(self.receptor_chains) == 1:
            content = content.replace("_receptorChainsNames_", 
                                      json.dumps(list(self.receptor_chains))) # represent as json list one string
        else: # receptor has 2 chains
            content = content.replace("_receptorChainsNames_", json.dumps(self.receptor_chains))
        
        # Remove leading spaces and write content into file
        content = textwrap.dedent(content)
        output_path = os.path.join(self.base_dir_control, "peEq.conf")
        with open(output_path, 'w') as file:
            file.write(content)


    def create_pele_conf_production(self):
        
        # Content for control file
        content = os.path.join(os.path.dirname(__file__), "pele_templates/peProd.conf")
        with open(content, "r") as file:
            content = file.read()
        
        # Get interacting residues
        # only interface between receptor chain 1 and ligand chain in case the receptor has 2 chains
        if len(self.receptor_chains) == 1 and self.receptor_chain_1 == None and self.receptor_chain_2 == None:
                analyzer = InterfaceAnalyzer(
                    pdb_path=self.pdb_path,
                    receptor_chain=self.receptor_chains,
                    ligand_chain=self.ligand_chain,
                )
        else: # receptor has 2 chains
            analyzer = InterfaceAnalyzer(
                pdb_path=self.pdb_path,
                receptor_chain=self.receptor_chain_1,
                ligand_chain=self.ligand_chain,
                )
                
        _ = analyzer.get_interaction_matrix()
        _ , _ = analyzer.expand_interface(neighborhood=3)
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
        content = content.replace("_ligandChainName_", self.ligand_chain)
        content = content.replace("_RIR_", json.dumps(interaction_matrix_receptor_chain_expanded_lst_pele))
        content = content.replace("_LIR_", json.dumps(interaction_matrix_ligand_chain_expanded_lst_pele))
        content = content.replace("_PDB_", self.pdb_file)
        
        if len(self.receptor_chains) == 1:
            content = content.replace("_receptorChainsNames_", 
                                      json.dumps(list(self.receptor_chains))) # represent as json list one string
        else: # receptor has 2 chains
            content = content.replace("_receptorChainsNames_", json.dumps(self.receptor_chains))
        
        # Remove leading spaces and write content into file
        content = textwrap.dedent(content)
        output_path = os.path.join(self.base_dir_control, "peProd.conf")
        with open(output_path, 'w') as file:
            file.write(content)