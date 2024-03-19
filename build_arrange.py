# Author: Curwen Pei Hong Tan
# Usage: python build_arrange.py

import os
ef_dataset_path = "/work/ratul1/chuen/Lanmodulin-data/ef_pdbs" # <-- change path here

files_in_ef = os.listdir(ef_dataset_path)
if ".DS_Store" in files_in_ef: files_in_ef.remove(".DS_Store")
files_in_ef = [f"{ef_dataset_path}/{x}" for x in files_in_ef]

motif_dict = {}

for file in files_in_ef:
    with open(file, 'r') as f:
        first_line = f.readline().strip()
        motif_type = first_line.split(":")[-1].strip()

        if motif_type in motif_dict:
            motif_dict[motif_type].append(os.path.basename(file))
        else:
            motif_dict[motif_type] = [os.path.basename(file)]

for motif_type, pdb_files in motif_dict.items():
    print(f"Motif Type: {motif_type}")
    print("PDB files: ")
    for pdb_file in pdb_files:
        print(pdb_file)
    print()


for key in motif_dict:
    print(key, len(motif_dict[key]))