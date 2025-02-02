import os 
from src.utils.logger_factory import LoggerFactory
from src.utils.os_utils import display_full_dataframe
import os
import sys
from pyrosetta import *
from pyrosetta.toolbox import *
import textwrap


logger = LoggerFactory.get_logger(__name__, "INFO")
display_full_dataframe()


def fast_relax_protocol(
    pdb_file, 
    scorefxn, 
    relax_protocol
    ):
    
    """
    Run and definea FastRelax protocol to a given pose.
    
    Parameters:
    pdb_file : str
        Path to the PDB file to relax.
    scorefxn : pyrosetta.ScoreFunction
        Score function to evaluate the pose.
    relax_protocol : pyrosetta.protocols.relax.FastRelax
        FastRelax protocol to apply to the pose.
    
    Returns:
    e_0 : float
        Initial energy of the pose.
    e_relaxed : float
        Energy of the pose after relaxation.
    """
    
    # Initial pose
    pose = pose_from_pdb(pdb_file)
    e_0 = scorefxn(pose)
    
    # Side-chain or side-chain and backbone relaxation
    #movemap = pyrosetta.rosetta.core.kinematics.MoveMap() # Remove if you want to relax the backbone
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
    """
    Run FastRelax protocol to a given PDB file.
    
    Parameters:
    pdb_file : str
        Path to the PDB file to relax.
        
    Returns:
        relax_report.tsv : file with the energies of the relaxed poses.
    """
    
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
    path_to_fastrelax_script,
    pdbs_dir
):
    """
    Generate a bash script to run FastRelax protocol in parallel over multiple PDB files.
    
    Parameters:
    path_to_fastrelax_script : str
        Path to the script that contains the FastRelax protocol.
    pdbs_dir : str
        Path to the directory with PDB files to relax.
    
    Returns:
        fastrelax.sh : bash script to run FastRelax protocol in parallel.
    
    """
    
    # Get PDBs to relax
    items = os.listdir(pdbs_dir)
    pdbs = [pdb for pdb in items if pdb.endswith('.pdb')]
    
    # Write string variable with greasy relaxjob info
    relax_job_greasy = "relaxjob.txt\n"+"\n".join([f"python {path_to_fastrelax_script} {os.path.join(os.path.abspath(pdbs_dir), pdb)}" for pdb in pdbs])
    relax_job_greasy = textwrap.indent(relax_job_greasy, '    ')
    
    content = f"""
    #!/bin/bash
    #SBATCH --job-name=__RJN__
    #SBATCH --output=logs/%x_%j.out
    #SBATCH --error=logs/%x_%j.err 
    #SBATCH --ntasks=__NPDBS__
    #SBATCH --cpus-per-task=1 
    #SBATCH --time=1:00:00
    #SBATCH --account=bsc72 
    #SBATCH --qos=gp_bscls

    module load greasy
    ml anaconda
    source activate /gpfs/scratch/bsc72/ismael/conda_envs/diffpele

    cat <<EOL > __COMMANDS__
    EOL

    greasy relaxjob.txt
    rm relaxjob.txt    
    """
    
    # Adjustments
    content = content.replace('__NPDBS__', str(len(pdbs)))
    content = content.replace('__COMMANDS__', relax_job_greasy)
    content = content.replace('__RJN__', os.path.basename(pdbs_dir)[:2]+"_rel")
        
    # Remove indentation
    content = textwrap.dedent(content)
    content = content.lstrip("\n")
    
    # Write content
    output_path = os.path.join(os.path.dirname(pdbs_dir), 'fastrelax.sh')
    with open(output_path, "w") as file:
        file.write(content)