import os, sys
from Bio.PDB import PDBParser, DSSP

def parse_pdb(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    return structure

def get_secondary_structure(structure):
    model = structure[0]  
    dssp = DSSP(model, pdb_file)
    
    # Create a dictionary to map DSSP output to secondary structure characters
    sec_struc_map = {
        'H': 'H',  # Alpha helix
        'B': 'E',  # Isolated beta-bridge
        'E': 'E',  # Extended strand
        'G': 'H',  # 3-helix (310 helix)
        'I': 'H',  # 5-helix (pi helix)
        'T': 'C',  # Turn
        'S': 'C',  # Bend
        '-': 'C'   # Coil
    }

    sec_structure_seq = ""
    for key in dssp.keys():
        residue_info = dssp[key]
        sec_structure = residue_info[2]
        sec_structure_seq += sec_struc_map.get(sec_structure, 'C')  

    return sec_structure_seq

def main(pdb_file):
    structure = parse_pdb(pdb_file)
    sec_structure_seq = get_secondary_structure(structure)
    print(f">{pdb_file}")
    print(sec_structure_seq)

pdb_file = sys.argv[1]
main(pdb_file)

