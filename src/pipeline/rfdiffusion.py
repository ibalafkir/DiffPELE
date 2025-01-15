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


class RFdiffContigs:
    
    def __init__(
        self,
        pdb_path: str,
        receptor_chains: str,
        ligand_chain: str,
        distance_cutOff: float
    ):
        
        """
        Initialize the RFdiffContigs class.
        This class pretends to obtain contigs from a system
        
        Parameters:
        pdb_path (str): Path to the PDB file.
        receptor_chains (str): Chain ID for the receptor protein chains.
                               Accepts one chain: 'A'
                               or 2 comma-sepparated: 'B'
        ligand_chain (str): Chain ID for the ligand protein chain.
                            Accepts one chain only
        distance_cutOff (float): Distance cutoff for the interface analysis
        """
        
        # Set attributes

        self.pdb_path = pdb_path
        self.distance_cutOff = distance_cutOff
        self.ligand_chain = ligand_chain
        
        if len(receptor_chains) == 1:
            self.receptor_chains = receptor_chains
        if len(receptor_chains) == 3:
            self.receptor_chains = receptor_chains.split(',')
            self.receptor_chain_1 = self.receptor_chains[0]
            self.receptor_chain_2 = self.receptor_chains[1]




    def _validate_inputs(
        self,
        pdb_path: str,
        receptor_chains: str,
        ligand_chain: str,
        distance_cutOff: float
    ):
        """
        Validate inputs for the RFdiffContigs class.
        """

        # Validate string inputs
        for var_name, var_value in [
            ("pdb_path", pdb_path),
            ("receptor_chain", receptor_chains),
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

        """
        # Validate chains in PDB
        atom_df = PdbDf(pdb_path)
        atom_df.get_atoms()
        atom_df_chains = atom_df.get_chains_id()
        if receptor_chains not in atom_df_chains:
            raise ValueError(f"{receptor_chains} is not a valid chain ID in {pdb_path}")
        if ligand_chain not in atom_df_chains:
            raise ValueError(f"{ligand_chain} is not a valid chain ID in {pdb_path}")
        if receptor_chains == ligand_chain:
            raise ValueError(f"Chain IDs should be different")
        if receptor_chains == ' ' or ligand_chain == ' ':
            raise ValueError(f"Chain IDs should not be empty")
        if len(receptor_chains) > 1 or len(ligand_chain) > 1:
            raise ValueError(f"Chain IDs should be single characters each")     
        """
        
        
    #@staticmethod
    def _get_selected_res(self):
        
        """
        Get the selected residues from the interface analysis.
        
        Returns:
        interaction_rec_lst (list): List of interacting residues in the receptor
                                    e.g. ['A1', 'A2', 'A3']
        interaction_lig_lst (list): List of interacting residues in the ligand
                                    e.g. ['B1', 'B2', 'B3']
                                         ['C1', 'C2', 'C3']
        """
        
        # One receptor chain - One ligand chain
        if len(self.receptor_chains) == 1:
            
            analyzer = InterfaceAnalyzer(
                pdb_path=self.pdb_path,
                receptor_chain=self.receptor_chains,
                ligand_chain=self.ligand_chain,
                distance_cutOff=self.distance_cutOff
            )
        
            _ = analyzer.get_interaction_matrix()
            interaction_rec_df , interaction_lig_df = analyzer.expand_interface(neighborhood=3)
            self.interaction_rec_lst, self.interaction_lig_lst = analyzer.get_interacting_residues(mode='e')
                    
        # Two receptor chains - One ligand chain
        if len(self.receptor_chains) == 2:
                
            # Analyze both interfaces per sepparate
            analyzer_1 = InterfaceAnalyzer(
                    pdb_path=self.pdb_path,
                    receptor_chain=self.receptor_chain_1,
                    ligand_chain=self.ligand_chain,
                    distance_cutOff=self.distance_cutOff
                )

            analyzer_2= InterfaceAnalyzer(
                    pdb_path=self.pdb_path,
                    receptor_chain=self.receptor_chain_2,
                    ligand_chain=self.ligand_chain,
                    distance_cutOff=self.distance_cutOff
                )

            _ = analyzer_1.get_interaction_matrix()
            interaction_rec_1_df , interaction_lig_df = analyzer_1.expand_interface(neighborhood=3)
            interaction_rec_1_lst, interaction_lig_lst_1 = analyzer_1.get_interacting_residues(mode='e')
            
            _ = analyzer_2.get_interaction_matrix()
            interaction_rec_2_df , interaction_lig_df = analyzer_2.expand_interface(neighborhood=3)
            interaction_rec_2_lst, interaction_lig_lst_2 = analyzer_2.get_interacting_residues(mode='e')
            
            self.interaction_rec_lst = interaction_rec_1_lst + interaction_rec_2_lst # sum the lists of interacting residues in receptor
                
            # Combine the lists of interacting residues in ligand
            interaction_lig_lst = interaction_lig_lst_1 + interaction_lig_lst_2
            interaction_lig_lst_uniq = list(set(interaction_lig_lst))
            self.interaction_lig_lst =  sorted(interaction_lig_lst_uniq, key=lambda x: int(x[1:]))
        
        return self.interaction_rec_lst, self.interaction_lig_lst
    
    def _get_and_filter_chunks(self, res_lst):
        
        """
        Get the chunks of residues from the selected residues.
        
        Parameters:
        res_lst (list): List of selected residues.
                        e.g. ['A1', 'A2', 'A3']
        """
        
        res_lst_nochainID = [int(i[1:]) for i in res_lst]
        
        result = []
        subgroup = []

        for res_num in res_lst_nochainID:
            if not subgroup or res_num == subgroup[-1] + 1:
                subgroup.append(res_num)
            else:
                result.append(subgroup)
                subgroup = [res_num]
                
        if subgroup:
            result.append(subgroup)
            
        # Deletes isolate residue res_numbers like 5 in [1,2,3], [5], [8, 9, 10] / THIS FITS BETTER IN chunk_filter (it works here though)
        for i in result:
            if len(i)==1:
                result.remove(i)
        
        # Filter chunks shorter than 3
        result_filter = [i for i in result if len(i) >= 3]
        
        return result_filter

    def _conver_list_to_contig(self, pdb_path, interaction_resnum_lst):
        """
        Gets a list of interacting residues after filtering and converts it
        
        Parameters:
        chain_id (str): Chain ID of the chain to analyze.
        pdb_path (str): Path to PDB file
        """
        
        # chain_id
        # chain_start
        # chain_end
        # chain_allresnum
        #

    
    
if __name__ == '__main__':
        
    contigs = RFdiffContigs(
        pdb_path='/home/ibalakfi/Desktop/testdiffpele/4POU_b.pdb',
        receptor_chains='B',
        ligand_chain='A',
        distance_cutOff=12.0
        )
    
    interaction_rec_lst, interaction_lig_lst = contigs._get_selected_res()
    interaction_rec_lst_rec_chunks = contigs._get_and_filter_chunks(interaction_rec_lst)
    interaction_lig_lst_lig_chunks = contigs._get_and_filter_chunks(interaction_lig_lst)
    print(interaction_rec_lst)
    print(interaction_lig_lst)
    print(interaction_rec_lst_rec_chunks)
    print(interaction_lig_lst_lig_chunks)
