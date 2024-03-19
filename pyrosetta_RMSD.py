
# for score function
from pyrosetta import *
init()
from pyrosetta.rosetta.core.scoring import CA_rmsd
from pyrosetta.toolbox import *
from pyrosetta.teaching import *

def getPairwiseRMSD(pdb1, pdb2):
    # Load the two protein complexes (poses) from PDB files
    pose1 = pose_from_pdb(pdb1)
    pose2 = pose_from_pdb(pdb2)
    # Calculate the RMSD between the two protein complexes
    rmsd_value = CA_rmsd(pose1, pose2)

    return rmsd_value

pdb1 = "ef_pdbs/ABS68055-1_unrelaxed_rank_005_alphafold2_ptm_model_5_seed_000_EF1.pdb"
pdb2 = "ef_pdbs/ACB32191-1_unrelaxed_rank_005_alphafold2_ptm_model_4_seed_000_EF1.pdb"

print(getPairwiseRMSD(pdb1, pdb2))