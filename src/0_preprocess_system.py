from utils.pdb_utils import PdbDf, PdbHandler, validate_and_get_chains
from utils.logger_factory import LoggerFactory
from utils.os_utils import display_full_dataframe
import argparse
import shutil
import os


logger = LoggerFactory.get_logger(__name__, "INFO")
display_full_dataframe()


def parse_args():
    
    script_desc= \
        "Prepare system for the DiffPELE pipeline."
    pdb_path_h = \
        "Path to the PDB file. \
        The PDB file must not contain resID gaps"
    receptor_chains_h = \
        "Chains of the receptor to be used in the pipeline. \
        Accepts up to two chains ID (e.g.: H and L antibody \
        chains). If inputting two chains, separate them with a comma."
    ligand_chain_h = \
        "Chains of the ligand to be used in the pipeline. \
        Accepts only one chain ID"
    
    parser = argparse.ArgumentParser(description = script_desc)
    
    parser.add_argument("--pdb_path",
                        "-pdb", 
                        type = str, 
                        required = True,
                        help = pdb_path_h)
    parser.add_argument('--receptor_chains',
                        '-rc', 
                        type=str, 
                        help=receptor_chains_h)
    parser.add_argument('--ligand_chain', 
                        '-lc',
                        type=str, 
                        help=ligand_chain_h)
    
    return parser.parse_args()
    
    
def main(pdb_path, receptor_chains, ligand_chain):

    # Parse PDB file and validate
    PdbHandler(pdb_path)
    logger.info(f"Parsing completed for the pipeline preprocessor: {os.path.realpath(pdb_path)}")

    # Validate and get chains
    chains_to_keep = validate_and_get_chains(receptor_chains, ligand_chain)
    logger.info(f"Chains to keep: {chains_to_keep}")
    
    # Create preprocessing directory and copy input pdb file
    pdb_path = os.path.realpath(pdb_path) # Make sure the path is absolute
    pdb_path_filename = os.path.basename(pdb_path) # E.g.: 1a2k.pdb
    pdb_path_name_noext = pdb_path_filename[:-4] # E.g.: 1a2k
    pdb_path_dir = os.path.dirname(pdb_path) # E.g.: /path/to/dir where dir holds 1a2k.pdb
    
    preprocessing_dir = os.path.join(pdb_path_dir, f"{pdb_path_name_noext}_preprocessed")
    os.makedirs(preprocessing_dir, exist_ok=True)
    
    pdb_path_new = os.path.join(preprocessing_dir, pdb_path_filename)
    shutil.copy(pdb_path, pdb_path_new)
    
    #### Main: saving chains to keep ####
    
    ## Define names for intermediate/temporal files
    
    pdb_step_1 = pdb_path_new[:-4] + "_noHETATM.pdb"
    pdb_step_2 = pdb_step_1[:-4] + "_selectedChains.pdb"
    pdb_step_3 = pdb_step_2[:-4] + "_tidied.pdb"
    pdb_step_4 = pdb_step_3[:-4] + "_insertionCodesRemoved.pdb"
    pdb_step_5 = pdb_step_4[:-4] + "_correctedTERs.pdb"
    pdb_step_6 = pdb_step_5[:-4] + "_tidied.pdb"
    
    ## Main
    PdbHandler.remove_hetatm(pdb_path_new)
    PdbHandler.save_chains(pdb_step_1, chains_to_keep)
    PdbHandler.tidy(pdb_step_2)
    PdbHandler.fix_insertion_codes(pdb_step_3)
    PdbHandler.fix_TERs(pdb_step_4)
    PdbHandler.tidy(pdb_step_5)
    
    ## Remove intermediate temporal files
    os.remove(pdb_path_new)
    os.remove(pdb_step_1)
    os.remove(pdb_step_2)
    os.remove(pdb_step_3)
    os.remove(pdb_step_4)
    os.remove(pdb_step_5)
    
    ## Rename final file
    os.rename(pdb_step_6, pdb_path_new[:-4]+"_preprocessed.pdb")    
    logger.info(f"Preprocessing completed for the pipeline preprocessor in: {os.path.realpath(preprocessing_dir)}")
    
if __name__ == "__main__":
    
    args = parse_args()
    main(
        args.pdb_path, 
        args.receptor_chains, 
        args.ligand_chain
        )