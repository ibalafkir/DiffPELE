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


class RFdiffusionContigs:
    
    def __init__(
        self,
        pdb_path: str,
        receptor_chains: str,
        ligand_chain: str,
        distance_cutOff: float = 12.0
    ):
        
        """
        Initialize the RFdiffusionContigs class.
        
        Parameters:
        pdb_path (str): Path to the PDB file.
        receptor_chains (str): Chain ID for the receptor protein chains.
        ligand_chain (str): Chain ID for the ligand protein chain.
        distance_cutOff (float): Distance cutoff for the RFdiffusion analysis.
        """
        
        # Set attributes
        self.pdb_path = pdb_path
        self.receptor_chains = receptor_chains
        self.ligand_chain = ligand_chain
        self.distance_cutOff = distance_cutOff
        
        self.chains_to_keep, self.receptor_chains, self.ligand_chain = validate_and_get_chains(
            self.receptor_chains, self.ligand_chain
        )
        
        # Define attributes for use in other methods
        
        
        
        # Preprocessing if necessary
    
    
    def _validate_inputs(
        
    ):
    
    
    #@staticmethod
    def _get_resnum(self, df):
        pass
    
    
    #@staticmethod
    def _get_chunks(self, resnum_lst):
        pass
    
    
    #@staticmethod
    def _filter_chunks(self, resnum_lst):
        pass
    
    
    #@staticmethod
    def get_contigs(self):
        """
        Gets the contigs between the interface of chain-chain (so they need to be specified).
        """
        # TODO figure out inputs in this method
        pass