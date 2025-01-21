from utils.pdb_utils import PdbDf, PdbHandler
from utils.logger_factory import LoggerFactory
from utils.os_utils import display_full_dataframe
from pipeline.pele import PeleSetup
import argparse
import shutil
import os


logger = LoggerFactory.get_logger(__name__, "INFO")
display_full_dataframe()


def parse_args():
    
    script_desc= \
        "Generate PELE production and equilibration protocols for single \
        or diffusion mode."
    pdb_path_h = \
        "Path to the PDB file. Necessary to compute the com_distance metric \
        of PELE. If using multiple mode, specify one of the multiple systems. If \
        running PELE with diffusion models, specify the original system"
    receptor_chains_h = \
        "Chains of the receptor. Accepts up to two \
        chains ID (e.g.: H and L antibody chains). \
        Input 2 chains by comma-sepparation."
    ligand_chain_h = \
        "Chains of the ligand. Accepts only one"
    distance_cutOff_h = \
        "Distance cut-off between the ligand and the receptor. \
        Default value is 12.0 Angstrom."
    mode_h = \
        "PELE simulation mode. Options: single, multiple (or diffusion, here it is the same). \
        Single mode: only one initial system. \
        Multiple mode: multiple initial systems. Used for diffusion \
        models"
    multiple_pdb_path_h = \
        "Path to the folder containing multiple PDB files. \
        Used for multiple mode. Not required for single mode."
    
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
    parser.add_argument('--distance_cutOff',
                        '-dc',
                        type=float, 
                        default=12.0,
                        help=distance_cutOff_h)
    parser.add_argument('--mode', 
                        '-m',
                        type=str,
                        default='single',
                        help=mode_h)
    parser.add_argument('--multiple_pdb_path',
                        '-mpdb',
                        type=str,
                        default = None,
                        help=multiple_pdb_path_h)
    
    return parser.parse_args()


def main(
    pdb_path: str,
    receptor_chains: str,
    ligand_chain: str,
    distance_cutOff: float,
    mode: str,
    multiple_pdb_path: str
    ):
    
    simulation_sys = PeleSetup(
        pdb_path = pdb_path,
        receptor_chains = receptor_chains,
        ligand_chain = ligand_chain,
        distance_cutOff = distance_cutOff
    )
    
    pdb_code = os.path.splitext(os.path.basename(pdb_path))[0]

    if mode == 'single':
        
        # Validate mode
        if multiple_pdb_path is not None:
            raise ValueError(
                "Multiple PDB path is not required for single mode."
            )
        
        # NBDsuite
        simulation_sys.create_nbdsuite_runner(
            suiteJobRunnerName=pdb_code[:2]+'_sui',
            nCPUs=1
        )
        simulation_sys.create_nbdsuite_input(
            nCPUs=1
        )
        
        # PELE equilibration
        simulation_sys.create_equilibration_runner(
            nCPUs=16,
            equiJobRunnerName=pdb_code[:2]+'_equi'
        )
        simulation_sys.create_control_adaptive_equilibration(
            outputPathName="outEQ",
            nEpochs=1,
            nSteps=40,
            nCPUs=16
        )        
        simulation_sys.create_pele_conf_equilibration()        

        # PELE production
        simulation_sys.create_production_runner(
            nCPUs=64,
            prodJobRunnerName=pdb_code[:2]+'_prod'
        )
        simulation_sys.create_control_adaptive_production(
            outputPathName="outPROD",
            nEpochs=1,
            nSteps=300,
            nCPUs=64
        )
        simulation_sys.create_pele_conf_production()
        
    if mode == 'multiple' or mode == 'diffusion' or mode == 'dif':
        
        # Valuidate mode
        if multiple_pdb_path is None:
            raise ValueError(
                "Specify the path to the folder containing \
                the multiple PDB initial systems."
            )       
            
        # Get the number of initial systems and the number of equilibrated systems
        # The production is run with the top 33% (total energy) of the equilibrated systems
        initial_pdbs = os.listdir(multiple_pdb_path)
        n_initial_pdbs = len(initial_pdbs) 
        n_equi_to_production = round(n_initial_pdbs*0.33)           
        
        # NBDsuite
        # Same number of CPUs as initial systems/PDBs
        simulation_sys.create_nbdsuite_runner(
            suiteJobRunnerName=pdb_code[:2]+'_d_sui',
            nCPUs = n_initial_pdbs 
        )
        simulation_sys.create_nbdsuite_input(
            nCPUs = n_initial_pdbs 
        )
        
        # PELE equilibration
        # 5 CPUs per initial system/PDB plus 1 CPU
        simulation_sys.create_equilibration_runner(
            nCPUs = 5*n_initial_pdbs + 1,
            equiJobRunnerName=pdb_code[:2]+'_d_eq'
        )
        simulation_sys.create_control_adaptive_equilibration(
            outputPathName="outEQ",
            nEpochs=1,
            nSteps=40,
            nCPUs=5*n_initial_pdbs + 1
        )        
        simulation_sys.create_pele_conf_equilibration()        

        # PELE production
        # 16 CPUs per initial system/PDB plus 1 CPU
        simulation_sys.create_production_runner(
            nCPUs = 16*n_equi_to_production + 1,
            prodJobRunnerName=pdb_code[:2]+'_d_pro'
        )
        simulation_sys.create_control_adaptive_production(
            outputPathName="outPROD",
            nEpochs=1,
            nSteps=200,
            nCPUs = 16*n_equi_to_production + 1,
        )
        simulation_sys.create_pele_conf_production()
        
        # -----------------------------------------------------------------
        # Once PELE files are generated, several adjustments need to be done
        # -----------------------------------------------------------------
        
        # Remove the initial systems/PDB in _pele/pdbs and copy those of multiple_pdb_path
        os.remove(
            os.path.join(
                os.path.abspath(pdb_path)[:-4]+'_pele', 'pdbs', pdb_code+'.pdb'
            )
        )
        
        # Copy the initial systems/PDB to _pele/pdbs
        for pdb in initial_pdbs:
            shutil.copy(
                os.path.join(multiple_pdb_path, pdb),
                os.path.join(
                    os.path.abspath(pdb_path)[:-4]+'_pele', 'pdbs', pdb
                )
            )
            
        # Change the PELE directory name to specify that it runs at multiple mode
        os.rename(
            os.path.abspath(pdb_path)[:-4]+'_pele',
            os.path.abspath(pdb_path)[:-4]+'_dif_pele'
        )
        logger.info("You have generated PELE simulation files for multiple mode (i.e. starting at multiple initial systems).")
        logger.info("Numbe of initial systems: "+str(n_initial_pdbs)+"; Number of equilibrated systems to run in production: "+str(n_equi_to_production))
        logger.info("Since this pipeline only starts at multiple initial systems using diffusion models, the PELE root directory name has been changed to specify that for making it clear.")
    
if __name__ == "__main__":
    
    args = parse_args()
    main(
        pdb_path = args.pdb_path,
        receptor_chains = args.receptor_chains,
        ligand_chain = args.ligand_chain,
        distance_cutOff = args.distance_cutOff,
        mode = args.mode,
        multiple_pdb_path = args.multiple_pdb_path
    )
