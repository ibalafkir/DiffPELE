import argparse
import os


def parse_args():
    
    script_desc = \
        "Rename atom names in a PDB file to match the REF15 force field."
    pdb_h = \
        "Path to the PDB file"
    
    parser = argparse.ArgumentParser(description='PELE poses processor')
    parser.add_argument('-pdb', type=str, help='Path to the PDB file')
    return parser.parse_args()


def postprocess(input_file):
    """
    Given a PDB, substitute atom names to match the REF15 force field
    
    Parameters:
    input_file (str): Path to the PDB file
    
    Returns:
    The processed PDB in its path with the same name but p- before it indicating 
    it has been processed
    """
    
    input_file_dir = os.path.dirname(input_file)
    input_file_basename = os.path.basename(input_file)
    out_file_name = os.path.join(input_file_dir,'p_'+input_file_basename)
    
    
    with open(input_file, 'r') as infile, open(out_file_name, 'w') as outfile:
        for line in infile:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                residue_name = line[17:20]
                if residue_name.strip() == "CYT" or residue_name.strip() == "CYX":
                    line = line[:17] + "CYS" + line[20:]
                if residue_name.strip() == "HIE" or residue_name.strip() == "HIP" or residue_name.strip() == "HID":
                    line = line[:17] + "HIS" + line[20:]      
            outfile.write(line)


def main(input_file):

    # Input is a single PDB
    if input_file.endswith('pdb'):

        postprocess(input_file)

    # Input is a dir
    else:
        
        items = os.listdir(input_file)
        pdbs = [pdb for pdb in items if pdb.endswith('.pdb')]

        for pdb in pdbs:
            postprocess(os.path.abspath(pdb))

if __name__ == "__main__":

    args = parse_args()
    pdb = args.pdb
    
    main(pdb)
    
    