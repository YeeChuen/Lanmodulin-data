# Author: Curwen Pei Hong Tan
# Usage: python build_arrange.py

import os
import numpy as np
import matplotlib.pyplot as plt

user_path = "C:/Users/ty_ch/Documents/Projects/ChowdhuryVSC/Lanmodulin-data" # <-- TODO: change to your path

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
print("example rmsd_dict value\n", rmsd_dict[('WP_124161803-1_unrelaxed_rank_005_alphafold2_ptm_model_3_seed_000_EF3.pdb', 'WP_225712291-1_unrelaxed_rank_005_alphafold2_ptm_model_5_seed_000_EF2.pdb')])
# NOTE: use the rmsd_dict for rmsd, you dont have to build your own dictionary. 

ef_dataset_path = f"{user_path}/ef_pdbs"

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

# NOTE 3/28/2024: need to produce heatmap per motif type
# Plot heatmap
def heatmapPlot(data, save_name, title = "halo"): # # <-- NOTE: having default here is okay, but you need to use it in plt.title
    plt.figure(figsize=(10, 8))
    #plt.imshow(data, cmap='viridis', origin='lower', interpolation='nearest')
    plt.imshow(data, cmap='viridis') # <-- NOTE: use this instead.
    plt.colorbar(label='RMSD')
    #plt.title(f'RMSD Heatmap {"halo"}') # TODO: implement so that the title shows motif type, replace "?" with your answer.
    plt.title(f'RMSD Heatmap {title}') # NOTE: this is the correct implementation.
    plt.xlabel('EF hand PDB File')
    plt.ylabel('EF hand PDB File')

    # NOTE: you don't have to know this, this just labels the heatmap ticks
    plt

    plt.savefig(save_name)

for motif_type in motif_dict:
    pdb_name_list = motif_dict[motif_type]
    data = np.zeros((len(pdb_name_list), len(pdb_name_list)))
    
    # Nested loop starts here
    
    # This is just an example, it is incorrect.
    for i, pdb1 in enumerate(files_in_ef):
        for j, pdb2 in enumerate(files_in_ef):
            pdb1_name = os.path.basename(pdb1)
            pdb2_name = os.path.basename(pdb2)
            '''
            if pdb1_name != pdb2_name:
                rmsd_value = rmsd_dict.get((pdb1_name, pdb2_name), None)
                if rmsd_value is not None:
                    heatmap_matrix[i, j] = rmsd_value
                else:
                    heatmap_matrix[i, j] = np.nan
            '''
    
    heatmapPlot(data, "example_heatmap.png", "Example Title")

    data = np.zeros((len(pdb_name_list), len(pdb_name_list)))

    # TODO: implement nested loop here
    for i, pdb1_name in enumerate(pdb_name_list):
        for j, pdb2_name in enumerate(pdb_name_list):
            if j <= i: continue
            value = 0
            try:
                value = rmsd_dict[(pdb1_name, pdb2_name)]
            except:
                value = rmsd_dict[(pdb2_name, pdb1_name)]

            data[i, j] = value
            data[j, i] = value

            '''# NOTE: below is not necessary actually 
            if i != j:
                rmsd_value = rmsd_dict.get((pdb1_name, pdb2_name), None)
                if rmsd_value is not None: 
                #if rmsd_value: # <-- NOTE: this is the way to do in python to check if a variable is not None.
                    data[i, j] = rmsd_value
                else: 
                    data[i, j] = np.nan
            '''

    print(motif_type)
    print(pdb_name_list)
    print(data)

    #heatmapPlot(data, f"heatmap_{motif_type}")
    heatmapPlot(data, f"heatmap_{motif_type}", title = motif_type) # <-- NOTE: you need to add title here also.