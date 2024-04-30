# Author: Curwen Pei Hong Tan

import os

ef_dataset_path = "C:/Users/ty_ch/Documents/Projects/ChowdhuryVSC/Lanmodulin-data/ef_pdbs"

files_in_ef = os.listdir(ef_dataset_path)
if ".DS_Store" in files_in_ef: files_in_ef.remove(".DS_Store")
files_in_ef = [f"{ef_dataset_path}/{x}" for x in files_in_ef]

motif_dict = {}

for file in files_in_ef:
    with open(file, 'r') as f:
        """
        third_line = f.readline().strip()
        motif_type = third_line.split(":")[-1].strip()
        """
        motif_type = None
        for i, line in enumerate(f):
            if i == 2:
                motif_type = line.strip().split(":")[-1].strip()
                break
        
        if motif_type:
            if motif_type in motif_dict:
                motif_dict[motif_type].append(os.path.basename(file))
            else:
                motif_dict[motif_type] = [os.path.basename(file)]
        else:
            print(f"Motif type not found: {file}")

from tqdm import tqdm

def getDictFromFile(file):
    sol = {}
    with open(file, 'r') as f:
        contents = f.readlines()
        contents = [x.strip() for x in contents if x.strip()] 

        print(f"Loading from {file}...")
        for x in tqdm(contents):
            # Splitting the line into tuple (pairs) and RMSD score
            pair, rmsd = x.split(':')
            pair = eval(pair)  # Convert the string representation of tuple to actual tuple
            rmsd = float(rmsd)

            sol[pair] = rmsd
            
    return sol

# Using the function to get the dictionary
rmsd_dict = getDictFromFile("ef_rmsd_result_v2.txt")
# NOTE: use the rmsd_dict for rmsd, you dont have to build your own dictionary. 

import os
import numpy as np
import matplotlib.pyplot as plt

def calculate_rmsd(pdb1, pdb2):
    # Placeholder function to calculate RMSD value between two PDB files
    # You need to implement this function
    pass

ef_dataset_path = "C:/Users/ty_ch/Documents/Projects/ChowdhuryVSC/Lanmodulin-data/ef_pdbs"

files_in_ef = os.listdir(ef_dataset_path)
if ".DS_Store" in files_in_ef: files_in_ef.remove(".DS_Store")
files_in_ef = [f"{ef_dataset_path}/{x}" for x in files_in_ef]

motif_dict = {}

for file in files_in_ef:
    with open(file, 'r') as f:
        motif_type = None
        for i, line in enumerate(f):
            if i == 2:
                motif_type = line.strip().split(":")[-1].strip()
                break
        
        if motif_type:
            if motif_type in motif_dict:
                motif_dict[motif_type].append(os.path.basename(file))
            else:
                motif_dict[motif_type] = [os.path.basename(file)]
        else:
            print(f"Motif type not found: {file}")

# Create a dictionary to store RMSD values
rmsd_values = {}

# Iterate over each motif type in the motif dictionary
for motif_type, pdb_files in motif_dict.items():
    for pdb1 in pdb_files:
        for pdb2 in pdb_files:
            pdb1_name = os.path.basename(pdb1)
            pdb2_name = os.path.basename(pdb2)
            if pdb1_name != pdb2_name:
                # Calculate RMSD value between pdb1 and pdb2
                rmsd_value = calculate_rmsd(pdb1, pdb2) # <-- NOTE: calculate_rmsd is not implemented
                rmsd_values[(pdb1_name, pdb2_name)] = rmsd_value
# NOTE: your rmsd_values dictionary is empty
print(rmsd_value)  # <-- see the result

# Create a matrix of RMSD values for heatmap
num_pdb_files = len(motif_dict[pdb_files[0]]) # <-- NOTE: This is wrong, motif_dict uses motif type as key, not pdb file name as key.
# NOTE: KeyError here. num_pdb_files is nothing.

heatmap_matrix = np.zeros((num_pdb_files, num_pdb_files))

for i, pdb1 in enumerate(pdb_files):
    for j, pdb2 in enumerate(pdb_files):
        pdb1_name = os.path.basename(pdb1)
        pdb2_name = os.path.basename(pdb2)
        if pdb1_name != pdb2_name:
            heatmap_matrix[i, j] = rmsd_values[(pdb1_name, pdb2_name)]

print(heatmap_matrix)

# Plot heatmap
plt.figure(figsize=(10, 8))
plt.imshow(heatmap_matrix, cmap='viridis', origin='lower', interpolation='nearest')
plt.colorbar(label='RMSD')
plt.title('RMSD Heatmap')
plt.xlabel('PDB File Index')
plt.ylabel('PDB File Index')
plt.show()