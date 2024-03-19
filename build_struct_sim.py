# Author: Yee Chuen Teoh
# python build_struct_sim.py

#______
# imports
import os
import numpy as np

# TMalign
from tmtools.io import get_structure, get_residue_data
from tmtools.testing import get_pdb_path
from tmtools import tm_align
from tmtools import TMScore

#______
# functions
# check for unique seq
def checkUnique():
    with open('EFhands.fasta', 'r') as f:
        contents = f.readlines()
        seq = set()
        for line in contents:
            if line[0] == ">":
                continue
            line = line.replace("\n", "")
            seq.add(line)
        
        print(len(seq))

def structureSimilarity():
    ef_dir_path = "C:/Users/ty_ch/Documents/Projects/ChowdhuryVSC/Lanmodulin-data/ef_pdbs"
    ef_path_list = os.listdir(ef_dir_path)
    ef_path_list = [f"{ef_dir_path}/{x}" for x in ef_path_list]

    s = get_structure(get_pdb_path(ef_path_list[0][:-4]))
    print(ef_path_list[0][:-4])
    chain = next(s.get_chains())
    coords1, seq1 = get_residue_data(chain)
    print(seq1)
    print(coords1.shape)

    s = get_structure(get_pdb_path(ef_path_list[8][:-4]))
    print(ef_path_list[8][:-4])
    chain = next(s.get_chains())
    coords2, seq2 = get_residue_data(chain)
    print(seq2)
    print(coords2.shape)

    result = tm_align(coords1, coords2, seq1, seq2)

    print(result.t)
    print(result.u)
    print(result.tm_norm_chain2)
    print(result.tm_norm_chain1)

    tmtools.TMScore(get_pdb_path(ef_path_list[0][:-4]), get_pdb_path(ef_path_list[8][:-4]))

#______
# main
def main():
    print(f"Unique sequence: {checkUnique()}")
    structureSimilarity()
    pass

if __name__ == "__main__":
    main()
