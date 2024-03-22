# Author: Yee Chuen Teoh
# Usage: python pyrosetta_RMSD.py

#______
# import
import os

# for score function
from pyrosetta import *
init()
from pyrosetta.rosetta.core.scoring import CA_rmsd
from pyrosetta.toolbox import *
from pyrosetta.teaching import *

#______
# function
def getPairwiseRMSD(pdb1, pdb2):
    # Load the two protein complexes (poses) from PDB files
    pose1 = pose_from_pdb(pdb1)
    pose2 = pose_from_pdb(pdb2)
    # Calculate the RMSD between the two protein complexes
    rmsd_value = CA_rmsd(pose1, pose2)

    return rmsd_value

def equalSequence(pdb1, pdb2):
    # Load the two protein complexes (poses) from PDB files
    pose1 = pose_from_pdb(pdb1)
    pose2 = pose_from_pdb(pdb2)

    return pose1.sequence() == pose2.sequence()

#______
# main
def main():
    ef_pdbs_path = "/work/ratul1/chuen/Lanmodulin-data/ef_pdbs"

    save_file = "/work/ratul1/chuen/Lanmodulin-data/ef_rmsd_result.txt"
    seq_gr2_file = "/work/ratul1/chuen/Lanmodulin-data/ef_equal_seq_rmsd_gr2.txt"
    seq_le2_file = "/work/ratul1/chuen/Lanmodulin-data/ef_equal_seq_rmsd_le2.txt"

    with open(save_file, 'w') as f:
        print(f"Created file: {save_file}")
    with open(seq_gr2_file, 'w') as f:
        print(f"Created file: {seq_gr2_file}")
    with open(seq_le2_file, 'w') as f:
        print(f"Created file: {seq_le2_file}")

    ef_pdbs_list = os.listdir(ef_pdbs_path)

    pose_dict = {}
    for x in ef_pdbs_list:
        pose_dict[f"{ef_pdbs_path}/{x}"] = pose_from_pdb(f"{ef_pdbs_path}/{x}")

    ef_pdbs_list = [f"{ef_pdbs_path}/{x}" for x in ef_pdbs_list]

    ef_dict = {}
    gr2_dict = {}
    le2_dict = {}

    for i, key1 in enumerate(list(pose_dict.keys())):
        print(f"currently on iter: {i}.")
        for j, key2 in enumerate(list(pose_dict.keys())):
            if j <= i: continue
            pose1 = pose_dict[key1]
            pose2 = pose_dict[key2]

            # v-- get the pdb file name rather than the whole path.
            i_name = key1.split("/")[-1]
            j_name = key2.split("/")[-1]
            
            # v-- basic checker to make sure this pairwise is unique. (SHOULD BE UNIQUE.)
            if (i_name, j_name) in ef_dict: raise ValueError("\nEF pair already added to dictionary.\nDuplicate should not happen.")
            
            rmsd = CA_rmsd(pose1, pose2) # <-- get rmsd score from pyRosetta
            if pose1.sequence() == pose2.sequence(): # <-- check if both pdb has the same sequence.
                if rmsd > 0.2:
                    gr2_dict[(i_name, j_name)] = rmsd # <-- this is what we're interested in, same seq, different structure.
                elif rmsd <= 0.2: 
                    le2_dict[(i_name, j_name)] = rmsd # <-- also save the other in case we needed it.
            
            ef_dict[(i_name, j_name)] = rmsd # <-- save all pairwise comparison RMSD value.

    with open(save_file, 'w') as f:
        '''
        for key in ef_dict:
            f.write(f"{key}:{ef_dict[key]}")
        '''
        f.write(str(ef_dict))
    with open(seq_gr2_file, 'w') as f:
        f.write(str(gr2_dict))
    with open(seq_le2_file, 'w') as f:
        f.write(str(le2_dict))

    with open(save_file.replace(".txt", "_v2.txt"), 'w') as f:
        for key in ef_dict:
            f.write(f"{key}:{ef_dict[key]}\n")
    with open(seq_gr2_file.replace(".txt", "_v2.txt"), 'w') as f:
        for key in gr2_dict:
            f.write(f"{key}:{gr2_dict[key]}\n")
    with open(seq_le2_file.replace(".txt", "_v2.txt"), 'w') as f:
        for key in le2_dict:
            f.write(f"{key}:{le2_dict[key]}\n")

    #pdb1 = "ef_pdbs/ABS68055-1_unrelaxed_rank_005_alphafold2_ptm_model_5_seed_000_EF1.pdb"
    #pdb2 = "ef_pdbs/ACB32191-1_unrelaxed_rank_005_alphafold2_ptm_model_4_seed_000_EF1.pdb"
    #print(getPairwiseRMSD(pdb1, pdb2))

if __name__ == "__main__":
    main()