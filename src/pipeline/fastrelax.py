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


import os
import sys
from pyrosetta import *
from pyrosetta.toolbox import *

def fast_relax(
    pdb_file, 
    scorefxn, 
    relax_protocol
    ):
    
    # Initial pose
    pose = pose_from_pdb(pdb_file)
    e_0 = scorefxn(pose)
    
    # Side-chain or side-chain and backbone relaxation
    movemap = pyrosetta.rosetta.core.kinematics.MoveMap()
    movemap.set_bb(False)  # do not allow backbone movements. Remove if you want to relax the backbone
    movemap.set_chi(True)  # allow side chain movements. Remove if you want to relax the backbone
    movemap.set_jump(False)  # disable rigid-body movements. Remove if you want to relax the backbone
    relax_protocol.set_movemap(movemap) # Remove if you want to relax the backbone
    relax_protocol.apply(pose)
    
    # Relax pose
    e_relaxed = scorefxn(pose)
    
    # File outputs
    logger.info(f'System {os.path.abspath(pdb_file)} relaxed from {e_0} to {e_relaxed}')
    pose.dump_pdb(f'{pdb_file[:-4]}_relax.pdb')
    
    return e_0, e_relaxed


if __name__ == '__main__':

    init('-relax:default_repeats 2')
    pdb_file = sys.argv[1]

    # Set scorer and FastRelax protocol
    max_relax_iter = 500
    scorefxn = get_fa_scorefxn()
    relax = pyrosetta.rosetta.protocols.relax.FastRelax()
    relax.set_scorefxn(scorefxn)
    relax.max_iter(max_relax_iter)

    # Apply FastRelax protocol
    logger.info("Relaxing:", pdb_file)
    e_0, e_relaxed = fast_relax(pdb_file, scorefxn, relax)
    
    # Write energies to a file
    tsv_string = 'pdb\tE_initial\tdelta_E\tE_final\n'
    tsv_string += f'{os.path.abspath(pdb_file)}\t{e_0}\t{e_relaxed-e_0}\t{e_relaxed}\n'
    with open(os.path.join(os.path.dirname(pdb_file), 'relax_report.tsv'), 'a') as f:
        f.write(tsv_string)        