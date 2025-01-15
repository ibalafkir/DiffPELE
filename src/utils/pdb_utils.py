import numpy as np
from Bio.PDB import PDBParser, Superimposer
from biopandas.pdb import PandasPdb
import copy
import pandas as pd
from pdbtools import pdb_fixinsert, pdb_tidy, pdb_selatom, pdb_sort, pdb_reatom, pdb_delchain, pdb_delelem, pdb_reres
import os


def validate_and_get_chains(receptor_chains, ligand_chain):

    """
    Validate and get the chains to keep in the PDB file.

    Parameters:
    receptor_chains (str): Chain ID for the receptor protein chains.
                           Accepts one chain: 'A'
                           or 2 comma-sepparated: 'B'
    ligand_chain (str): Chain ID for the ligand protein chain.
                        Accepts one chain only

    Returns:
    chains_to_keep (list): List of chains to keep in the PDB file
    receptor_chains (str or list): String for one chain, list for two
    ligand_chain (str): Ligand chain
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
    else: # len(receptor_chains) == 3, meaning
        receptor_chains = receptor_chains.split(',')
        chains_to_keep.append(receptor_chains[0])
        chains_to_keep.append(receptor_chains[1])
        chains_to_keep.append(ligand_chain)

    return chains_to_keep, receptor_chains, ligand_chain


class PdbDf:
    """
    Manipulates PDB files using Pandas DataFrames.
    """
    
    def __init__(
        self, 
        pdb_path: str
    ):
        """
        Initializes the PdbDf object by loading the PDB file.
        
        Parameters:
        pdb_path (str): Path to the PDB file.
        """
        self.pdb = pdb_path
        self.ppdb_atoms = None
        self.ppdb_chains_id = None
        self.ppdb_atoms_ca = None
        
        
        self._validate_inputs(
            pdb_path = self.pdb
        )
        self.ppdb = PandasPdb().read_pdb(self.pdb)
        
        
    def _validate_inputs(
        self,
        pdb_path: str
    ):
        """
        Validates the input parameters.
        Only validations regarding PDB input are performed
        since chain and residue validations are done in scripts
        """
        
        if not isinstance(pdb_path, str):
            raise TypeError(f"Expected string, got {type(pdb_path)}")
        
        if not os.path.exists(pdb_path):
            raise FileNotFoundError(f"File not found: {pdb_path}")
 

    def get_atoms(self):
        """
        Extracts ATOM information in a pandas DataFrame.
        
        Returns:
        pd.DataFrame: DataFrame containing ATOM records.
        """
        self.ppdb_atoms = self.ppdb.df['ATOM']
        return self.ppdb_atoms

    
    def get_chains_id(self):
        """
        Extracts unique chain IDs from the ATOM records.
        
        Returns:
        list: List of unique chain IDs.
        """
        self.ppdb_chains_id = self.ppdb_atoms['chain_id'].unique().tolist()
        return self.ppdb_chains_id 

       
    def get_atoms_ca(self):
        """
        Extracts CA atoms in a pandas DataFrame.
        
        Returns:
        pd.DataFrame: DataFrame containing CA atoms.
        """
        self.ppdb_atoms_ca = self.ppdb_atoms[self.ppdb_atoms['atom_name'] == 'CA']
        return self.ppdb_atoms_ca


    @staticmethod
    def get_atoms_chain(
        atom_df: pd.DataFrame, 
        chain_id: str
    ):
        """
        Extracts atoms from a specific chain.
        
        Parameters:
        atom_df (pd.DataFrame): DataFrame containing ATOM records.
        chain_id (str): Chain ID.
        
        Returns:
        pd.DataFrame: DataFrame containing ATOM records from the specified chain.
        """
        
        def _validate_inputs():
            
            if not isinstance(atom_df, pd.DataFrame):
                raise TypeError(f"Expected pd.DataFrame, got {type(atom_df)}")
            
            if not isinstance(chain_id, str):
                raise TypeError(f"Expected string, got {type(chain_id)}")
            
            if chain_id not in atom_df['chain_id'].unique():
                raise ValueError(f"Chain ID not found: {chain_id}")
            
        _validate_inputs()
        atom_df_chain_id = atom_df[atom_df['chain_id'] == chain_id]
        
        return atom_df_chain_id


    @staticmethod
    def get_atom_relevant(
        atom_df: pd.DataFrame
    ):
        """
        Extracts relevant positional data from the ATOM records.
        
        Returns:
        pd.DataFrame: DataFrame containing positional data.
        """

        coordinates_info = [
            'chain_id', 
            'residue_number', 
            'residue_name', 
            'x_coord', 
            'y_coord', 
            'z_coord'
        ]
        
        atom_coordinates = atom_df[coordinates_info]
        
        return atom_coordinates
    

class PdbHandler:
    """
    General PDB processing functions that generate new PDBs 
    by applying these functions.
    Only validations regarding PDB input are performed
    since chain and residue validations are done in scripts
    """

    def __init__(
        self, 
        pdb_path: str
    ):
        """
        Initializes the PdbHandler object by loading the PDB file.
        
        Parameters:
        pdb_path (str): Path to the PDB file.
        """
        self.pdb = pdb_path
        
        self._validate_inputs(
            pdb_path = self.pdb
        )


    def _validate_inputs(
        self,
        pdb_path: str
    ):
        """
        Validates the input parameters.
        """
        
        if not isinstance(pdb_path, str):
            raise TypeError(f"Expected string, got {type(pdb_path)}")
        
        if not os.path.exists(pdb_path):
            raise FileNotFoundError(f"File not found: {pdb_path}")


    @staticmethod
    def fix_insertion_codes(
        pdb_path: str
    ):
        """
        Fixes insertion codes in the PDB file.
        
        Parameters:
        pdb_path (str): Path to the PDB file.
        
        Returns:
        None: This method creates a new PDB file with the fixed content.
        """
        out_path = pdb_path[:-4]+'_insertionCodesRemoved.pdb' # Output name
        
        f = open(pdb_path, "rt") # Read original PDB
        f_new = open(out_path, "wt") # Write fixed PDB
        
        lines = f.readlines() # Read lines
        
        # Call the fixinsert function
        for modified_line in pdb_fixinsert.run(lines, []):
            f_new.write(modified_line)
        
        # Close files
        f.close()
        f_new.close()
    
    
    @staticmethod
    def fix_TERs(
        pdb_path: str
        ):
        """
        Fixes TER mistakes in the PDB file.

        Parameters:
        pdb_path (str): Path to the PDB file.
        
        Returns:
        None: This method creates a new PDB file with the fixed content.
        """
        out_path = pdb_path[:-4]+'_correctedTERs.pdb'
        
        with open(pdb_path, "r") as f_input, open(out_path, "w") as f_output:
            
            # Corrects TER mistakes
            for line in f_input:
                if line.startswith("TER"): 
                    f_output.write("TER\n")
                    atom_index = line.find("ATOM")
                    if atom_index != -1:
                        f_output.write(line[atom_index:])  
                else:
                    f_output.write(line)
    
    
    @staticmethod              
    def remove_hetatm(pdb_path):
        """
        Removes HETATM lines and their ANISOU lines (for now, other removals might be added here)
    
        Parameters:
        pdb_path (str): Path to the PDB file.
        
        Returns:
        None: This method creates a new PDB file with the new content.
        """
        with open(pdb_path, 'r') as f_input, open(pdb_path[:-4]+'_noHETATM.pdb', "w") as f_output:
            c = 0 # Counter
            for line in f_input: # Start reading
                if line.startswith('HETATM'): 
                    c += 1 # If bumping into a HETATM, counter stops being 0 -> this serves to know we've entered the HETATM area
                    continue # We don't write HETATM lines
                if c != 0 and line.startswith('ANISOU'): # If we are in the HETATM area (c = 0) and we bump into ANISOU lines
                    continue # We don't write them
                f_output.write(line)


    @staticmethod
    def remove_elements(pdb_path, elements_lst):
        """
        Removes all lines starting by all the elements indicated in the list 
        (e.g. ['ATOM', 'HETATM']) in the output PDB
        
        Parameters:
        pdb_path (str): Path to the PDB file.
        lst (list): List of elements to remove.
        
        Returns:
        None: This method creates a new PDB file with the non-removed elements.
        """
        with open(pdb_path, 'r') as f_input, open(pdb_path[:-4]+'_removed.pdb', "w") as f_output:
            for line in f_input:
                if not any(line.startswith(element) for element in elements_lst):
                    f_output.write(line)
    
    
    @staticmethod
    def save_chains(pdb_path, chains_lst):
        """
        Keeps the atom names from the desired chain names
        
        Parameters:
        pdb_path (str): Path to the PDB file.
        
        Returns:
        None: This method creates a new PDB file with the saved chains.
        """
        out_path = pdb_path[:-4]+'_selectedChains.pdb'
        
        pdb_df = PdbDf(pdb_path)
        atom_df = pdb_df.get_atoms()
        atom_df_ch = atom_df[atom_df['chain_id'].isin(chains_lst)]
        atom_df_ch['segment_id']=''
        
        # TODO Maybe change this
        new_pdb = PandasPdb().read_pdb(pdb_path)
        new_pdb.df['ATOM'] = atom_df_ch
        new_pdb.to_pdb(path=out_path, records=['ATOM'], gz = False)
    
    
    @staticmethod
    def renumber_atoms(pdb_path, first_atom_idx=1):
        """
        ---
        BUG TODO
        Troubles joining lines, test 4pou.pdb
        ---
        
        Renumbers a PDB starting from indicated atom number
        
        Parameters:
        pdb (str): Path to the PDB file.
        atom_num (int): Number to start renumbering from.
        
        Returns:
        None: This method creates a new PDB file with renumbered atoms.
        """
        out_path = pdb_path[:-4]+'_atomSorted.pdb'

        f = open(pdb_path, 'rt')
        f_new = open(out_path, 'wt')
        
        lines = f.readlines()        
        for modified_line in pdb_reatom.run(lines, first_atom_idx):
            f_new.write(modified_line)
            
        f.close()
        f_new.close()
    
    
    @staticmethod
    def renumber_residues(pdb_path, first_res_idx=1):
        """
        Renumbers a PDB starting from the desired residue number
        
        Warning: renumbers all residues regardless of the chain
                 and also numbers TER
        
        Parameters:
        pdb (str): Path to the PDB file.
        first_res_idx (int): Number to start renumbering from.
        
        Returns:
        None: This method creates a new PDB file with renumbered residues.
        """
        f = open(pdb_path, 'rt')
        f_resnum = open(pdb_path[:-4]+'_residueRenumbered.pdb', 'wt')
        
        lines = f.readlines()        
        for modified_line in pdb_reres.run(lines, first_res_idx):
            f_resnum.write(modified_line)
            
        f.close()
        f_resnum.close()


    @staticmethod
    def delete_atoms(pdb_path, atom_lst):
        """
        Deletes all atomic elements indicated in the list
        
        Warning: to delete C but not CA, input 'C ' (with the space)
        
        Parameters:
        pdb (str): Path to the PDB file.
        atom_lst (list): List of atomic elements to delete.
                         e.g.: ['C', 'O', 'N']
                         Some will depend on the FF used in the PDB
    
        Returns:
        None: This method creates a new PDB file with the new content.
        """
        out_path = pdb_path[:-4]+'_atomRemoved.pdb'
        
        f = open(pdb_path, 'rt')
        f_deleted = open(out_path, 'wt')
        
        lines = f.readlines()        
        for modified_line in pdb_delelem.run(lines, atom_lst):
            f_deleted.write(modified_line)
            
        f.close()
        f_deleted.close()


    @staticmethod
    def get_backbone(pdb_path):
        """
        Deletes all atoms but the ones belonging to the backbone
        Does not affect other lines (headers...)
        
        Parameters:
        pdb_path (str): Path to the PDB file.
        
        Returns:
        None: This method creates a new PDB file with the backbone residues.
        """
        
        out_path = pdb_path[:-4]+'_backbone.pdb'
        
        f = open(pdb_path, 'rt')
        f_backbone = open(out_path, 'wt')
        
        lines = f.readlines()        
        for modified_line in pdb_selatom.run(lines, [
            'CA','C','N','O'
            ]):
            f_backbone.write(modified_line)
            
        f.close()
        f_backbone.close()
        
   
    @staticmethod
    def get_atom(pdb_path):
        """
        Keeps only atom lines in the written PDB
        It does not write TER, END or any other type of line
        To solve this, use pdbtools-tidy function

        Parameters:
        pdb_path (str): Path to the PDB file.
        
        Returns:
        None: This method creates a new PDB file with atom lines only.
        """
        pdb_atom = pdb_path[:len(pdb_path)-4]+'_atom.pdb'
        ppdb = PandasPdb().read_pdb(pdb_path)
        ppdb.to_pdb(path= pdb_atom, records = ['ATOM'], gz=False)


    @staticmethod
    def tidy(pdb_path):
        """
        Detects chain ID and res ID changes to assign TER and END lines
        E.g.:
        it will add a TER between atom lines in which resid is x and the
        next x+y (y=>2)
        it will add a TER between atom lines where the chain ID is different 
        
        
        Parameters:
        pdb (str): Path to the PDB file.
        
        Returns:
        None: This method creates a new PDB file with the tidied content.
        """
        out_path = pdb_path[:-4]+'_tidied.pdb'
        
        f = open(pdb_path, 'rt')
        f_new = open(out_path, 'wt')
        
        lines = f.readlines()        
        for modified_line in pdb_tidy.run(lines, False):
            f_new.write(modified_line)
        
        f.close()
        f_new.close()


    @staticmethod
    def sort(pdb_path):
        """ 
        Sorts a PDB according to chain ID and residue number

        Parameters:
        pdb (str): Path to the PDB file.
        
        Returns:
        None: This method creates a new PDB file with the sorted content.
        """
        out_path = pdb_path[:-4]+'_sorted.pdb'
        
        f = open(pdb_path, 'rt')
        f_sorted = open(out_path, 'wt')
        
        lines = f.readlines()        
        for modified_line in pdb_sort.run(lines, 'CR'): # 'CR' orders first by chain and then by residues
            f_sorted.write(modified_line)
            
        f.close()
        f_sorted.close()
        
        