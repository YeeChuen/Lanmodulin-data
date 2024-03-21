# Author: Curwen Pei Hong Tan
# Usage: python report_ef_data.py

#______
# imports
import json

#______
# functions

#______
# main
def main():
    ef_rmsd_file = "ef_rmsd_result.txt"
    ef_gr2_rmsd_file = "ef_equal_seq_rmsd_gr2.txt"
    ef_le2_rmsd_file = "ef_equal_seq_rmsd_le2.txt"

    with open(ef_gr2_rmsd_file, 'r') as f:
        contents = f.readlines()
        print(len(contents))
        print("\n" in contents[0])

        if len(contents) != 1: raise ValueError("\nContent should be only 1.")
        if "\n" in contents[0]: raise ValueError("\nContent should not have new line.")

        rmsd_dict = eval(contents[0])

        print(len(rmsd_dict))
        print(rmsd_dict[('WP_124161803-1_unrelaxed_rank_005_alphafold2_ptm_model_3_seed_000_EF3.pdb', 'WP_246790065-1_unrelaxed_rank_005_alphafold2_ptm_model_5_seed_000_EF3.pdb')])
        


if __name__ == "__main__":
    main()
