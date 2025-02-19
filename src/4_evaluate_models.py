from utils.pdb_utils import PdbDf, PdbHandler
from utils.logger_factory import LoggerFactory
from utils.os_utils import display_full_dataframe
from pipeline.pele import PeleSetup
import argparse
import shutil
import os
from DockQ.DockQ import load_PDB, run_on_all_native_interfaces


"""
This script evaluates the quality of a model using DockQ.
For systems with 2 receptor chains, the chains must be merged into a single chain.
Strongly inspired by DockQ authors: https://github.com/bjornwallner/DockQ/issues/33
"""


logger = LoggerFactory.get_logger(__name__, "INFO")
display_full_dataframe()


def parse_args():
    
    script_desc= \
        "Run DockQ to assess the quality between a native PDB and a model. \
        The ligand must have 1 chain whereas the receptor can have 1 or 2."
    pdb_path_h = \
        "Path to the PDB file"
    model_path_h = \
        "Path to the model file"
    receptor_chains_h = \
        "Two chains ID comma-sepparated (e.g.: H,L) for the receptor protein or \
        one chain ID for the receptor protein."
    ligand_chain_h = \
        "Accepts only one chain ID"
    
    parser = argparse.ArgumentParser(description = script_desc)
    
    parser.add_argument("--pdb_path",
                        "-pdb", 
                        type = str, 
                        required = True,
                        help = pdb_path_h)
    parser.add_argument("--model_path",
                        "-mp", 
                        type = str, 
                        required = True,
                        help = model_path_h)
    parser.add_argument('--receptor_chains',
                        '-rc', 
                        type=str, 
                        help=receptor_chains_h)
    parser.add_argument('--ligand_chain', 
                        '-lc',
                        type=str, 
                        help=ligand_chain_h)
    
    return parser.parse_args()


def merge_chains(model, chains_to_merge):
    for chain in chains_to_merge[1:]:
        for i, res in enumerate(model[chain]):
            res.id = (chain, res.id[1], res.id[2])
            model[chains_to_merge[0]].add(res)
        model.detach_child(chain)
    model[chains_to_merge[0]].id = "".join(chains_to_merge)
    return model


def main(pdb_path, 
         model_path,
         receptor_chains, 
         ligand_chain):
    
    if len(receptor_chains) == 1:

        model = load_PDB(model_path)
        native = load_PDB(pdb_path)
        chain_map = {f"{receptor_chains}":f"{receptor_chains}", f"{ligand_chain}":f"{ligand_chain}"}
        dockq_dict = run_on_all_native_interfaces(model, native, chain_map=chain_map)
        first_key = next(iter(dockq_dict[0]))  
        values = dockq_dict[0][first_key]
        print("DockQ:", round(values['DockQ'],3), "iRMSD:", round(values['iRMSD'],3), "LRMSD:", round(values['LRMSD'],3), "fnat:", round(values['fnat'], 3), "fnonnat:", round(values['fnonnat'],3))

    if len(receptor_chains) == 3 and ',' in receptor_chains:
        
        model = load_PDB(model_path)
        native = load_PDB(pdb_path)
        receptor_chains = receptor_chains.split(',')
        model_merged = merge_chains(model, receptor_chains)
        native_merged = merge_chains(native, receptor_chains)
        chain_map = {f"{receptor_chains[0]}{receptor_chains[1]}":f"{receptor_chains[0]}{receptor_chains[1]}", f"{ligand_chain}":f"{ligand_chain}"}
        dockq_dict = run_on_all_native_interfaces(model_merged, native_merged, chain_map=chain_map)
        first_key = next(iter(dockq_dict[0]))
        values = dockq_dict[0][first_key]
        print("DockQ:", round(values['DockQ'],3), "iRMSD:", round(values['iRMSD'],3), "LRMSD:", round(values['LRMSD'],3), "fnat:", round(values['fnat'], 3), "fnonnat:", round(values['fnonnat'],3))


if __name__ == "__main__":
    
    args = parse_args()
    
    main(
        pdb_path = args.pdb_path,
        model_path = args.model_path,
        receptor_chains = args.receptor_chains,
        ligand_chain = args.ligand_chain
    )