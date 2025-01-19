import os 
from src.utils.logger_factory import LoggerFactory
from src.utils.os_utils import display_full_dataframe
import os
import sys
from pyrosetta import *
from pyrosetta.toolbox import *


logger = LoggerFactory.get_logger(__name__, "INFO")
display_full_dataframe()


def fast_relax_protocol(
    pdb_file, 
    scorefxn, 
    relax_protocol
    ):
    
    # Initial pose
    pose = pose_from_pdb(pdb_file)
    e_0 = scorefxn(pose)
    
    # Side-chain or side-chain and backbone relaxation
    movemap = pyrosetta.rosetta.core.kinematics.MoveMap() # Remove if you want to relax the backbone
    #movemap.set_bb(False)  # do not allow backbone movements. Remove if you want to relax the backbone
    #movemap.set_chi(True)  # allow side chain movements. Remove if you want to relax the backbone
    #movemap.set_jump(False)  # disable rigid-body movements. Remove if you want to relax the backbone
    #relax_protocol.set_movemap(movemap) # Remove if you want to relax the backbone
    relax_protocol.apply(pose)
    
    # Relax pose
    e_relaxed = scorefxn(pose)
    
    # File outputs
    logger.info(f'System {os.path.abspath(pdb_file)} relaxed from {e_0} to {e_relaxed}')
    pose.dump_pdb(f'{pdb_file[:-4]}_relax.pdb')
    
    return e_0, e_relaxed


def run_fast_relax(
    pdb_file
):
    
    # Initialize PyRosetta
    init('-relax:default_repeats 2')
    
    # Set scorer and FastRelax protocol
    max_relax_iter = 500
    scorefxn = get_fa_scorefxn()
    relax = pyrosetta.rosetta.protocols.relax.FastRelax()
    relax.set_scorefxn(scorefxn)
    relax.max_iter(max_relax_iter)
    
    # Apply FastRelax protocol
    logger.info("Relaxing: "+pdb_file)
    e_0, e_relaxed = fast_relax_protocol(pdb_file, scorefxn, relax)
    
    # Write energies to a file
    tsv_string = 'pdb\tE_initial\tdelta_E\tE_final\n'
    tsv_string += f'{os.path.abspath(pdb_file)}\t{e_0}\t{e_relaxed-e_0}\t{e_relaxed}\n'
    with open(os.path.join(os.path.dirname(pdb_file), 'relax_report.tsv'), 'a') as f:
        f.write(tsv_string)        

    return e_0, e_relaxed


def make_fastrelax_runner(
    fast_relax_script,
    pdbs_dir
):
    content = """
    #!/bin/bash
    #SBATCH --job-name=relax
    #SBATCH --output=logs/%x_%j.out
    #SBATCH --error=logs/%x_%j.err 
    #SBATCH --ntasks=__NPDBS__
    #SBATCH --cpus-per-task=1 
    #SBATCH --time=1:00:00
    #SBATCH --account=bsc72 
    #SBATCH --qos=gp_bscls


    module load greasy
    ml anaconda
    source activate /gpfs/projects/bsc72/conda_envs/pyrosetta

    cat <<EOL > relaxjob.txt
    python <path/to/fast_relax_script> <iterate_over_pdbs_in_pdb_dir>
    EOL

    greasy relaxjob.txt
    rm relaxjob.txt
    """
    
    
if __name__ == '__main__':

    # Parse arguments
    pdb_file = sys.argv[1]
    run_fast_relax(pdb_file)