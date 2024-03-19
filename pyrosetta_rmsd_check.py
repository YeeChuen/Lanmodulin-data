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

if __name__ == "__main__":
    main()