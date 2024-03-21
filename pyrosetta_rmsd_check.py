# Author: Yee Chuen Teoh
'''
Usage: 
    python pyrosetta_rmsd_check.py --pdb1='WP_124161803-1_unrelaxed_rank_005_alphafold2_ptm_model_3_seed_000_EF3.pdb' --pdb2='WP_225712291-1_unrelaxed_rank_005_alphafold2_ptm_model_5_seed_000_EF2.pdb'

('WP_225712291-1_unrelaxed_rank_005_alphafold2_ptm_model_5_seed_000_EF2.pdb', 'WP_231327943-1_unrelaxed_rank_005_alphafold2_ptm_model_2_seed_000_EF3.pdb'): 0.24888131022453308
    python pyrosetta_rmsd_check.py --pdb1='WP_225712291-1_unrelaxed_rank_005_alphafold2_ptm_model_5_seed_000_EF2.pdb' --pdb2='WP_231327943-1_unrelaxed_rank_005_alphafold2_ptm_model_2_seed_000_EF3.pdb'

'''

#______
# import
import os
import argparse

# for score function
from pyrosetta import *
init()
from pyrosetta.rosetta.core.scoring import CA_rmsd
from pyrosetta.toolbox import *
from pyrosetta.teaching import *

#______
# function
def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb1', type=str, required=True)
    parser.add_argument('--pdb2', type=str, required=True)

    args = parser.parse_args()
    return args

def getPairwiseRMSD(pdb1, pdb2):
    # Load the two protein complexes (poses) from PDB files
    pose1 = pose_from_pdb(pdb1)
    pose2 = pose_from_pdb(pdb2)
    # Calculate the RMSD between the two protein complexes
    rmsd_value = CA_rmsd(pose1, pose2)



    return rmsd_value

#______
# main
def main():
    args = parser()
    print(getPairwiseRMSD(f"./ef_pdbs/{args.pdb1}", f"./ef_pdbs/{args.pdb2}"))

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
        if i == 5: break
        
        print(f"{str(i)}/{str(len(list(pose_dict.keys())))}")
        print(pose_dict[key1].pdb_info())

        for j, key2 in enumerate(list(pose_dict.keys())):
            if j <= i: continue # <-- skip if j is less than or equal to i, since it will be repeated info.

            rmsd_value = CA_rmsd(pose_dict[key1], pose_dict[key2])

            print(key1)
            print(key2)
            print(rmsd_value)



if __name__ == "__main__":
    main()