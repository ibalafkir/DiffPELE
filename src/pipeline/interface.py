import numpy as np
import pandas as pd
import os 
import copy as cp
from src.utils.pdb_utils import PdbDf, PdbHandler
from src.utils.logger_factory import LoggerFactory
from src.utils.os_utils import display_full_dataframe
import argparse


logger = LoggerFactory.get_logger(__name__, "INFO")
display_full_dataframe()


class InterfaceAnalyzer:
    """
    Analyze protein-protein interfaces in a PDB file.
    
    The PDB file needs to have this format nuances:
    - No alt_locs (A, B, C...)
    - No jumps in residue numbers (e.g. Ser21 and after Lys23)
    
    This class uses methods that use the PdbDf class to extract
    atom information
    
    It can only analyze interfaces between two single chains
    """
    
    
    def __init__(
        self,
        pdb_path: str,
        receptor_chain: str,
        ligand_chain: str,
        distance_cutOff: float = 12.0
    ):
        """
        Initialize the InterfaceAnalyzer class.
        
        Parameters:
        pdb_path (str): Path to the PDB file.
        receptor_chain (str): Chain ID for the first protein chain.
        ligand_chain (str): Chain ID for the second protein chain.
        distance_cutOff (float): Distance cut-off for interaction analysis. Default is 12.0.
        """

        # Set attributes
        self.pdb_path = pdb_path
        self.receptor_chain = receptor_chain
        self.ligand_chain = ligand_chain
        self.distance_cutOff = distance_cutOff
        
        # Define attributes for later use
        self.interaction_matrix = None
        self.atom_df_CA_receptor_chain_re = None
        self.atom_df_CA_ligand_chain_re = None
        self.interaction_receptor_chain_df = None
        self.interaction_ligand_chain_df = None     
        
        # Validate inputs
        self._validate_inputs(
            pdb_path=pdb_path,
            receptor_chain=receptor_chain,
            ligand_chain=ligand_chain,
            distance_cutOff=distance_cutOff
        )   
        
        
    def _validate_inputs(
        self,
        pdb_path: str,
        receptor_chain: str,
        ligand_chain: str,
        distance_cutOff: float
    ):
        """
        Validate the inputs for the class.
        """

        # Validate string inputs
        for var_name, var_value in [
            ("pdb_path", pdb_path),
            ("receptor_chain", receptor_chain),
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
        if receptor_chain not in atom_df_chains:
            raise ValueError(f"{receptor_chain} is not a valid chain ID in {pdb_path}")
        if ligand_chain not in atom_df_chains:
            raise ValueError(f"{ligand_chain} is not a valid chain ID in {pdb_path}")
        if receptor_chain == ligand_chain:
            raise ValueError(f"Chain IDs should be different")
        if receptor_chain == ' ' or ligand_chain == ' ':
            raise ValueError(f"Chain IDs should not be empty")
        if len(receptor_chain) > 1 or len(ligand_chain) > 1:
            raise ValueError(f"Chain IDs should be single characters each")
    
    
    def get_interaction_matrix(self):
        
        """
        Get a DataFrame with the interacting residues between two chains
        in a PDB file. The interacting residues are defined as residues
        whose C-alpha atoms are closer than a distance cut-off.
        
        Returns:
        interaction_matrix (pd.DataFrame): The interacting residues between the two chains
        """
        
        # Read PDB as pandas DataFrame
        pdb = PdbDf(self.pdb_path) # validates PDB and opens/reads all lines
        atom_df = pdb.get_atoms() # Get atoms as DataFrame, deletes others
        
        # Get CA atoms
        atom_df_CA = pdb.get_atoms_ca() # deletes other atoms while keeping only CAs

        # Get atoms for chain 1 and chain 2
        atom_df_CA_receptor_chain = PdbDf.get_atoms_chain(atom_df_CA, self.receptor_chain)
        atom_df_CA_ligand_chain = PdbDf.get_atoms_chain(atom_df_CA, self.ligand_chain) 
        
        # Get relevant ('re', for shorter variable names) data from atoms. In each chain
        atom_df_CA_receptor_chain_re = PdbDf.get_atom_relevant(atom_df_CA_receptor_chain)
        atom_df_CA_ligand_chain_re = PdbDf.get_atom_relevant(atom_df_CA_ligand_chain)
        
        # Calculating distances between CAs
        atom_df_CA_receptor_chain_re_coords = atom_df_CA_receptor_chain_re[["x_coord", "y_coord", "z_coord"]].values
        atom_df_CA_ligand_chain_re_coords = atom_df_CA_ligand_chain_re[["x_coord", "y_coord", "z_coord"]].values
        distances = np.linalg.norm(
            atom_df_CA_receptor_chain_re_coords[:, None] - atom_df_CA_ligand_chain_re_coords,
            axis=2
        )
        distances_cutOff = np.where(distances <= self.distance_cutOff)
        atom_df_CA_receptor_chain_re_interactingRes = atom_df_CA_receptor_chain_re.iloc[distances_cutOff[0]]
        atom_df_CA_ligand_chain_re_interactingRes = atom_df_CA_ligand_chain_re.iloc[distances_cutOff[1]]       

        # Create distance matrix
        self.interaction_matrix = pd.DataFrame({
            f'Residues - Chain {self.receptor_chain}': atom_df_CA_receptor_chain_re_interactingRes['residue_name'].values,
            f'Residue number - Chain {self.receptor_chain}': atom_df_CA_receptor_chain_re_interactingRes['residue_number'].values,
            'X': atom_df_CA_receptor_chain_re_interactingRes['x_coord'].values,
            'Y': atom_df_CA_receptor_chain_re_interactingRes['y_coord'].values,
            'Z': atom_df_CA_receptor_chain_re_interactingRes['z_coord'].values,
            f'Residues - Chain {self.ligand_chain}': atom_df_CA_ligand_chain_re_interactingRes['residue_name'].values,
            f'Residue number - Chain {self.ligand_chain}': atom_df_CA_ligand_chain_re_interactingRes['residue_number'].values,
            'X': atom_df_CA_ligand_chain_re_interactingRes['x_coord'].values,
            'Y': atom_df_CA_ligand_chain_re_interactingRes['y_coord'].values,
            'Z': atom_df_CA_ligand_chain_re_interactingRes['z_coord'].values,
            'Distance between CAs': distances[distances_cutOff]
        })
        
        # Define self attributes for later use in next methods #
        # ---------------------------------------------------- #
        
        # Drop duplicates serves to delete alt_locs and only take into account A (e.g. AMET, not BMET, CMET...)
        self.atom_df_CA_receptor_chain_re = atom_df_CA_receptor_chain_re #.drop_duplicates()
        self.atom_df_CA_ligand_chain_re = atom_df_CA_ligand_chain_re #.drop_duplicates()
    
        # Drop duplicates and sort by residue number to get the interacting residues in order and without repetitions
        self.interaction_receptor_chain_df = atom_df_CA_receptor_chain_re_interactingRes.drop_duplicates().sort_values(by='residue_number')
        self.interaction_ligand_chain_df = atom_df_CA_ligand_chain_re_interactingRes.drop_duplicates().sort_values(by='residue_number')
            
        # Return interaction matrix
        return self.interaction_matrix
    
    
    def _expand_resnum_lst(self, resLst_toExpand, 
                           firstRes, lastRes, 
                           neighborhood):
        
        """
        Expand a list of residue numbers by a neighborhood.
        Private method used in expand_interface.
        """
        
        
        resLst_Expanded = cp.deepcopy(resLst_toExpand)
        
        # Iterate over resLst_toExpand and if there is a difference <= neighborhood 
        # between 2 following numbers, add the missing numbers in result
        for i in range(len(resLst_toExpand)-1):
            if resLst_toExpand[i+1] - resLst_toExpand[i] <= neighborhood:
                for j in range(resLst_toExpand[i]+1, resLst_toExpand[i+1]):
                    resLst_Expanded.append(j)
        
        # Sort the list and delete repeated
        resLst_Expanded = list(set(resLst_Expanded))
        resLst_Expanded = sorted(resLst_Expanded)
        
        # If first res of the list - firstres of the chain is <= neighborhood, add the missing numbers
        if resLst_Expanded[0] - firstRes <= neighborhood:
            for i in range(firstRes, resLst_Expanded[0]):
                resLst_Expanded.append(i)
        
        # Sort the list and delete repeated
        resLst_Expanded = list(set(resLst_Expanded))
        resLst_Expanded = sorted(resLst_Expanded)
        
        # If the lastres - last res of the lst of the chain is <= neighborhood, add the missing numbers
        if lastRes - resLst_Expanded[-1] <= neighborhood:
            for i in range(resLst_Expanded[-1]+1, lastRes+1):
                resLst_Expanded.append(i)
        
        # Sort the list and delete repeated
        resLst_Expanded = list(set(resLst_Expanded))
        resLst_Expanded = sorted(resLst_Expanded)
        
        return resLst_Expanded
        
        
    def expand_interface(self, neighborhood):
        
        """
        Expand the interacting residues by a neighborhood.
        
        Returns:
        interaction_receptor_chain_df_expanded (pd.DataFrame): The expanded interacting residues pandas dataframe for chain 1.
        interaction_ligand_chain_df_expanded (pd.DataFrame): The expanded interacting residues pandas dataframe for chain 2.
        """
        
        # Prepare residue data
        residues_receptor_chain = self.atom_df_CA_receptor_chain_re['residue_number'].values.tolist()
        residues_ligand_chain = self.atom_df_CA_ligand_chain_re['residue_number'].values.tolist()
        residues_receptor_chain_firstRes = residues_receptor_chain[0]
        residues_receptor_chain_lastRes = residues_receptor_chain[-1]
        residues_ligand_chain_firstRes = residues_ligand_chain[0]
        residues_ligand_chain_lastRes = residues_ligand_chain[-1]
        residues_receptor_chain_interactingRes = list(
            self.interaction_receptor_chain_df[
                    'residue_number'
                    ].values
                )
        residues_ligand_chain_interactingRes = list(
                    self.interaction_ligand_chain_df[
                        'residue_number'
                        ].values
                )

        # Expand interaction residue numbers lists by neighborhood
        residues_receptor_chain_interactingRes_expanded = self._expand_resnum_lst(
                                                       residues_receptor_chain_interactingRes, 
                                                       residues_receptor_chain_firstRes, 
                                                       residues_receptor_chain_lastRes, 3
                                                       )
        residues_ligand_chain_interactingRes_expanded = self._expand_resnum_lst(
                                                        residues_ligand_chain_interactingRes,
                                                        residues_ligand_chain_firstRes,
                                                        residues_ligand_chain_lastRes, 3
                                                        )

        # In the interaction dataframes, keep only the residue numbers that are in the expanded lists
        self.interaction_receptor_chain_df_expanded = self.atom_df_CA_receptor_chain_re[
            self.atom_df_CA_receptor_chain_re['residue_number'].isin(set(residues_receptor_chain_interactingRes_expanded))
        ]
        self.interaction_ligand_chain_df_expanded = self.atom_df_CA_ligand_chain_re[
            self.atom_df_CA_ligand_chain_re['residue_number'].isin(set(residues_ligand_chain_interactingRes_expanded))
        ]
        
        # Output the expanded by neighborhood interaction dataframes
        return self.interaction_receptor_chain_df_expanded, self.interaction_ligand_chain_df_expanded
        

    def get_interacting_residues(self, mode):
        
        """
        Get the interacting residues in a PDB file.
        
        Returns:
        receptor_chain_interaction_resnum_chainID (list): List of interacting residues for chain 1.
        ligand_chain_interaction_resnum_chainID (list): List of interacting residues for chain 2.
        """
        
        if mode in ['expanded', 'exp', 'e', 'E']:
            
            # Transform the expanded residue interaction dataframes 
            # to lists ['chainID'+'residue_number_1', ...]
            
            receptor_chain_interaction_chainID = list(
                                     (self.interaction_receptor_chain_df_expanded[
                                         'chain_id'
                                         ].values)
                                     )
            receptor_chain_interaction_resnum = list((self.interaction_receptor_chain_df_expanded[
                                        'residue_number'
                                        ].values)
                                    )      
            receptor_chain_interaction_resnum_chainID = [str(receptor_chain_interaction_chainID[i]) + 
                                          str(receptor_chain_interaction_resnum[i]) 
                                          for i in range(len(receptor_chain_interaction_resnum))
                                          ]
            ligand_chain_interaction_chainID = list(
                                        (self.interaction_ligand_chain_df_expanded[
                                            'chain_id'
                                            ].values)
                                        )
            ligand_chain_interaction_resnum = list(
                                        (self.interaction_ligand_chain_df_expanded[        
                                            'residue_number'
                                            ].values)
                                        )
            ligand_chain_interaction_resnum_chainID = [str(ligand_chain_interaction_chainID[i]) +
                                            str(ligand_chain_interaction_resnum[i])
                                            for i in range(len(ligand_chain_interaction_resnum))
                                            ]
        
        else:
            
            receptor_chain_interaction_chainID = list(
                                     (self.interaction_receptor_chain_df[
                                         'chain_id'
                                         ].values)
                                     )
            receptor_chain_interaction_resnum = list((self.interaction_receptor_chain_df[
                                        'residue_number'
                                        ].values)
                                    )      
            receptor_chain_interaction_resnum_chainID = [str(receptor_chain_interaction_chainID[i]) + 
                                          str(receptor_chain_interaction_resnum[i]) 
                                          for i in range(len(receptor_chain_interaction_resnum))
                                          ]
            ligand_chain_interaction_chainID = list(
                                        (self.interaction_ligand_chain_df[
                                            'chain_id'
                                            ].values)
                                        )
            ligand_chain_interaction_resnum = list(
                                        (self.interaction_ligand_chain_df[        
                                            'residue_number'
                                            ].values)
                                        )
            ligand_chain_interaction_resnum_chainID = [str(ligand_chain_interaction_chainID[i]) +
                                            str(ligand_chain_interaction_resnum[i])
                                            for i in range(len(ligand_chain_interaction_resnum))
                                            ]

        return receptor_chain_interaction_resnum_chainID, ligand_chain_interaction_resnum_chainID
        
        
############################################################################################################################################################

# Parsing args and main call for THIS script, which pretends to analyze the interface between 2 given chains of a system
    
def parse_args():
    
    script_desc= \
        "Analyze the interface between 2 chain of a given PDB."
    pdb_path_h = \
        "Path to the PDB file. \
        The PDB file must not contain resID gaps."
    receptor_chain_h = \
        "Chains of the receptor to be used in the pipeline. \
        Accepts only one chain ID (e.g.: H or L antibody \
        chain)."
    ligand_chain_h = \
        "Chains of the ligand to be used in the pipeline. \
        Accepts only one chain ID"
    distance_cutOff_h = \
        "Distance cut-off for interaction analysis. Default is 12.0. \
        The interacting residues are defined as residues whose C-alpha atoms \
        are closer than the distance cut-off."
    
    parser = argparse.ArgumentParser(description = script_desc)
    
    parser.add_argument("--pdb_path",
                        "-pdb", 
                        type = str, 
                        required = True,
                        help = pdb_path_h)
    parser.add_argument('--receptor_chain',
                        '-rc', 
                        type=str, 
                        help=receptor_chain_h)
    parser.add_argument('--ligand_chain', 
                        '-lc',
                        type=str, 
                        help=ligand_chain_h)
    parser.add_argument('--distance_cutOff',
                        '-dc',
                        type=float,
                        default=12.0,
                        help=distance_cutOff_h)
    
    return parser.parse_args()
    
    
def main(
    pdb_path: str,
    receptor_chain: str,
    ligand_chain: str,
    distance_cutOff: float = 12.0
    ):
    
    # Initialize the InterfaceAnalyzer class
    analyzer = InterfaceAnalyzer(
        pdb_path, 
        receptor_chain, 
        ligand_chain, 
        distance_cutOff
        )
    
    # Get the interaction matrix and write to csv
    interaction_matrix = analyzer.get_interaction_matrix()
    interaction_matrix.to_csv(f'{pdb_path[:-4]}+_interaction_matrix.csv', index=False)

    # Expand the interacting residues dataframes by a neighborhood and write to csv
    interaction_matrix_receptor_chain_expanded, interaction_matrix_ligand_chain_expanded \
        = analyzer.expand_interface(neighborhood=3)
    interaction_matrix_receptor_chain_expanded.to_csv(f'{pdb_path[:-4]}_interaction_matrix_expanded.csv', index=False)
    interaction_matrix_ligand_chain_expanded.to_csv(f'{pdb_path[:-4]}_interaction_matrix_expanded.csv.', header='False', mode='a', index=False)
    
    # Get the expanded interacting residues as lists
    interaction_matrix_receptor_chain_expanded_lst, interaction_matrix_ligand_chain_expanded_lst \
        = analyzer.get_interacting_residues(mode='expanded')
    
    # Include output pymol automatic selection
    interaction_matrix_receptor_chain_expanded_lst_noChainName = \
    [item.lstrip(receptor_chain) for item in interaction_matrix_receptor_chain_expanded_lst] # remove chain ID from list
    interaction_matrix_ligand_chain_expanded_lst_noChainName = \
    [item.lstrip(ligand_chain) for item in interaction_matrix_ligand_chain_expanded_lst] # remove chain ID from list
    interaction_matrix_receptor_chain_expanded_lst_noChainName_pymol = \
        "+".join(interaction_matrix_receptor_chain_expanded_lst_noChainName) # join list into string (1+2+3+...)
    interaction_matrix_ligand_chain_expanded_lst_noChainName_pymol = \
        "+".join(interaction_matrix_ligand_chain_expanded_lst_noChainName) # join list into string (1+2+3+...)
    pymol_selection = (
        f"select interface, (chain {receptor_chain} and resi {interaction_matrix_receptor_chain_expanded_lst_noChainName_pymol}) "
        f"or (chain {ligand_chain} and resi {interaction_matrix_ligand_chain_expanded_lst_noChainName_pymol})"
    )                                                                 # pymol selection string
    
    # Log the results
    logger.info("The interacting residues for chain 1 are: %s", interaction_matrix_receptor_chain_expanded_lst)
    logger.info("The interacting residues for chain 2 are: %s", interaction_matrix_ligand_chain_expanded_lst)
    logger.info("Run this command in the pymol command line to select residues from the interface: %s", pymol_selection)
    logger.info(f"The interaction matrix has been saved to {pdb_path[:-4]}_interaction_matrix.csv")
    logger.info(f"The expanded interaction matrix has been saved to {pdb_path[:-4]}_interaction_matrix_expanded.csv")
    
        
if __name__ == "__main__":
    
    args = parse_args()
    
    main(
        args.pdb_path, 
        args.receptor_chain, 
        args.ligand_chain, 
        args.distance_cutOff
        )
    
############################################################################################################################################################
