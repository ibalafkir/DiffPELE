from utils.pdb_utils import PdbDf, PdbHandler
from utils.logger_factory import LoggerFactory
from utils.os_utils import display_full_dataframe
from pipeline.fastrelax import fast_relax_protocol, run_fast_relax, make_fastrelax_runner
import argparse
import os


logger = LoggerFactory.get_logger(__name__, "INFO")
display_full_dataframe()


def parse_args():
    
    script_desc= \
        """
        Apply a FastRelax minimization protocol sequentally to all PDBs in a directory \
        or generate bash runners for running in parallel over multiple PDB files in the \
        same directory.
        
        Usage:
        
        To relax one PDB
        python ./src/2_fastrelax.py -m sequentally -pdb ./data/pdb_files/1a2k.pdb
        
        To relax multiple PDBs of a directory sequentally
        python ./src/2_fastrelax.py -m sequentally -pdb ./data/pdb_files/
        
        To generate a bash runner to run parallely over multiple PDBs in a directory
        python ./src/2_fastrelax.py -m parallel -pdb ./data/pdb_files/
        """
    mode = \
        "Mode of operation: 'sequentally' for running a single PDB file or 'parallel' \
        generating runner to apply parallely for multiple PDB files"
    pdb_path_h = \
        "Path to the directory with PDB files"

    parser = argparse.ArgumentParser(description = script_desc)
    
    parser.add_argument("--mode",
                        "-m", 
                        type = str, 
                        required = True,
                        choices = ['sequentally', 'parallel'],
                        help = mode)
    
    parser.add_argument("--pdb_path",
                        "-pdb", 
                        type = str, 
                        required = True,
                        help = pdb_path_h)
    
    return parser.parse_args()


def main():
    
    args = parse_args()
    
    if args.mode == 'sequentally':
        
        if args.pdb_path.endswith('.pdb'):
            run_fast_relax(args.pdb_path)
        else:
            items = os.listdir(args.pdb_path)
            pdbs = [pdb for pdb in items if pdb.endswith('.pdb')]

            for pdb in pdbs:
                run_fast_relax(os.path.join(args.pdb_path, pdb))
    
    elif args.mode == 'parallel':
            
        make_fastrelax_runner(
            path_to_fastrelax_script=f"{os.path.dirname(__file__)}/2_fastrelax.py -m sequentally -pdb",
            pdbs_dir=args.pdb_path
            )

if __name__ == "__main__":
    main()
