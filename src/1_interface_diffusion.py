from utils.pdb_utils import PdbDf, PdbHandler
from utils.logger_factory import LoggerFactory
from utils.os_utils import display_full_dataframe
from pipeline.pele import PeleSetup
from pipeline.rfdiffusion import RFdiffContigs, RFdiffFix, RFdiffSetUp
import argparse
import shutil
import os


logger = LoggerFactory.get_logger(__name__, "INFO")
display_full_dataframe()


def parse_args():
    
    script_desc= \
        "Generate RFdiffusion runner suitable for a cluster, in our case MNV."
    mode_h = \
        "This script can both generate the runner (and its inputs) but also correct \
        the output of the RFdiffusion. Options: setup, fix."
    pdb_path_h = \
        "Path to the PDB file. Necessary for both generation RFdiffusion runner and \
        diffusion models correction"
    receptor_chains_h = \
        "Chains of the receptor. Accepts up to two chains ID (e.g.: H and L antibody chains). \
        Input 2 chains by comma-sepparation. Not necessary if running the script in fix mode."
    ligand_chain_h = \
        "Chains of the ligand. Accepts only one. Not necessary if running the script in fix mode."
    distance_cutOff_h = \
        "Distance cut-off between the ligand and the receptor. \
        Default value is 12.0 Angstrom. Not necessary if running the script in fix mode."
    run_inference_path_h = \
        "Path to the run_inference RFdiffusion script. \
        Not necessary if running the script in fix mode."
    diffusion_models_input = \
        "Path to the diffusion models (path/to/diffusion_models/) OR \
        a single diffusion model (path/to/diffusion_model.pdb). This is \
        automatically handled by the class. Not necessary if running the script in setup mode."

    parser = argparse.ArgumentParser(description = script_desc)
    
    parser.add_argument("--mode",
                        "-m",
                        type = str,
                        required = True,
                        help = mode_h)
    parser.add_argument("--pdb_path",
                        "-pdb", 
                        type = str, 
                        required = True,
                        help = pdb_path_h)
    parser.add_argument('--receptor_chains',
                        '-rc', 
                        type=str, 
                        required=False,
                        default=None,
                        help=receptor_chains_h)
    parser.add_argument('--ligand_chain',
                        '-lc',
                        type=str, 
                        required=False,
                        default=None,
                        help=ligand_chain_h)
    parser.add_argument('--distance_cutOff',
                        '-dc',
                        type=float, 
                        required=False, 
                        default=12.0,
                        help=distance_cutOff_h)
    parser.add_argument('--run_inference_path',
                        '-rip',
                        type=str,
                        default='/gpfs/projects/bsc72/Repos/RFdiffusion/scripts/run_inference.py',
                        help=run_inference_path_h)
    
    parser.add_argument('--diffusion_models_input',
                        '-dmi',
                        type=str,
                        default=None,
                        help=diffusion_models_input)
    
    return parser.parse_args()


def main(
    mode: str,
    pdb_path: str,
    receptor_chains: str,
    ligand_chain: str,
    distance_cutOff: float,
    run_inference_path: str,
    diffusion_models_input: str
):
    
    if mode == 'setup':
        
        logger.info("Setting up the RFdiffusion runner.")
        
        # Small input validation
        if pdb_path is None:
            raise ValueError("PDB path is necessary for the setup mode.")
        if receptor_chains is None or ligand_chain is None:
            raise ValueError("Receptor and ligand chains are necessary for the setup mode.")
        if distance_cutOff is None:
            raise ValueError("Distance cut-off is necessary for the setup mode.")
        if run_inference_path is None:
            raise ValueError("Path to the run_inference script is necessary for the setup mode.")
        
        if diffusion_models_input is not None:
            logger.warning("Diffusion models input is not necessary for the setup mode.")
        
        # Organize directory
        path_to_pdb = os.path.abspath(pdb_path)        
        base_diff_dir = os.path.abspath(path_to_pdb)[:-4]+"_diffusion"
        os.makedirs(base_diff_dir)
        shutil.copy(pdb_path, base_diff_dir)
        pdb_path = os.path.join(base_diff_dir, os.path.basename(pdb_path))

        # Get contigs
        system = RFdiffContigs(
            pdb_path = pdb_path,
            receptor_chains = receptor_chains,
            ligand_chain = ligand_chain,
            distance_cutOff = distance_cutOff
        )
        system_contigs = system.get_contigs()
        with open(os.path.join(pdb_path[:-4]+"_contigs.txt"), "w") as f:
            if not system_contigs.endswith("\n"):
                f.write(system_contigs + "\n")
            else:
                f.write(system_contigs)
        
        # Make runner
        setup = RFdiffSetUp(
            pdb_path = pdb_path,
            run_inference_path = run_inference_path
            )
        rfdiffusion_runner = setup.make_runner()
        with open(os.path.join(base_diff_dir, "rfdiffusion_partial.sh"), "w") as f:
            f.write(rfdiffusion_runner)        

    elif mode == 'fix':

        logger.info("Fixing the RFdiffusion output.")
        
        # Small input validation
        if pdb_path is None:
            raise ValueError("PDB path is necessary for the fix mode.")
        if diffusion_models_input is None:
            raise ValueError("Diffusion models input is necessary for the fix mode.")
        for i in [
            receptor_chains,
            ligand_chain,
            distance_cutOff,
            run_inference_path
        ]:
            if i is not None:
                logger.warning(f"Input {i} is not necessary for the fix mode, but it does not affect correction")

        # Correct diffusion models, output is in the same directory as diffusion models
        diffusion_models_correction = RFdiffFix(
            original_pdb_path = pdb_path,
            diffusion_models_input = diffusion_models_input
        )
        diffusion_models_correction.correct_diff_models()
        
        # Make directory for the corrected diffusion models
        diffusion_models_input_path = os.path.dirname(diffusion_models_input)
        diffusion_models_input_name = os.path.basename(diffusion_models_input)
        fix_diffusion_models_path = os.path.join(diffusion_models_input_path, diffusion_models_input_name+"_fix")
        os.makedirs(fix_diffusion_models_path, exist_ok=True)
        
        # Move corrected PDBs to new directory
        items = os.listdir(diffusion_models_input)
        pdbs_to_move = [os.path.join(diffusion_models_input,pdb) for pdb in items if pdb.endswith('_fix.pdb')]
        for pdb in pdbs_to_move:
            shutil.move(pdb, fix_diffusion_models_path)
        
    else:
        raise ValueError("Mode not recognized. Please choose between 'setup' or 'fix'.")

if __name__ == "__main__":

    args = parse_args()
    main(
        mode = args.mode,
        pdb_path = args.pdb_path,
        receptor_chains = args.receptor_chains,
        ligand_chain = args.ligand_chain,
        distance_cutOff = args.distance_cutOff,
        run_inference_path = args.run_inference_path,
        diffusion_models_input = args.diffusion_models_input
    )


        

        