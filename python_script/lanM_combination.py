# Author: Yee Chuen Teoh
'''
usage:
    python lanM_combination.py --fasta "../fasta/Lanmodulins.fasta" --select "../ef_selected/stableEnergy_result.txt"

update:
creation of script (5/3/2024)
    - this script is created to combine the 9 fragment of selected fragments. 
'''
#_____
# imports
import argparse
import os
from tqdm import tqdm

# lanM_get_fragments import
from lanM_get_fragments import getFragments

#_____
# functions
def parser():
    parser = argparse.ArgumentParser(description="Read a PDF file from the RAG model.")
    parser.add_argument("--fasta", help=":The fasta file containing lanM, structure formated from the python script `build_ef_data.py`", required = True)
    parser.add_argument("--select", help=":A text file containing selected lanM and the part of EF obtain from the script `stableEnergy.py`", required = True)
    return parser.parse_args()

def parseSelectFile(select_file):
    with open(select_file, "r") as f:
        contents = f.readlines()
        contents = [x.replace("\n", "") for x in contents if x != "\n"]
        contents = [x.replace("\n", "") for x in contents if x != "\n"]
        contents = [("_".join(x.split("_")[:-10]), x.split("_")[-1].replace(".pdb", "")) for x in contents if x != "\n"]
        
        return contents

def makeHalfDesign(frags, parts, id):
    parts = [(tuple[0], tuple[1].lower()) for tuple in parts]
    half_design_seq = ""
    description = ""
    lanM_list = []
    ef_list = []
    for part in parts:
        if part[0] not in lanM_list:
            lanM_list.append(part[0])
        if "ef" in part[1]:
            ef_list.append(part[1])
        description += f"|{part}"
        half_design_seq += frags[part[0]][part[1]]

    description = ">" + str(id) + "_(" + ",".join(ef_list) + ")_(" + ",".join(lanM_list) + ")|" + description[1:]

    return (description, half_design_seq)

def buildParts(selects, build, start_i = 1, end_i = 3):
    print(f"\n--- Build parts with information: {build} ---")
    first_half = build[start_i]
    second_half = build[end_i]

    parts_list = []

    first_half_c = 0
    second_half_c = 0
    total_c = 0
    for tuple in selects:
        if tuple[1] == first_half:
            parts_1 = []
            for i in range(0, end_i):
                parts_1.append((tuple[0], build[i].strip().lower()))

            second_half_c = 0
            first_half_c += 1

            for tuple_j in selects:
                if tuple_j[1] == second_half:
                    parts_2 = []
                    for j in range(end_i, len(build)):
                        parts_2.append((tuple_j[0], build[j].strip().lower()))
                    second_half_c += 1
                    total_c += 1
                    
                    parts_list.append(parts_1 + parts_2)

    print(f"Number of first half (parts)({first_half}): {first_half_c}")
    print(f"Number of first half (parts)({second_half}): {second_half_c}")
    print(f"Number of combination (parts): {total_c}, {len(parts_list)}")

    return parts_list

def writeFastaTupleListToFile(fasta_list, save_name):
    with open(save_name, "w") as f:
        for tuple in fasta_list:
            f.write(f"{tuple[0]}\n")
            f.write(f"{tuple[1]}\n")

def combine(fasta_file, select_file):
    frags = getFragments(fasta_file)
    selects = parseSelectFile(select_file)

    test_half_design_parts = [('WP_238286702-1', 'head'),('WP_238286702-1', 'ef2'),('WP_238286702-1', 'linker2'), ('WP_074826049-1', 'ef3'), ('WP_074826049-1', 'tail')]
    makeHalfDesign(frags, test_half_design_parts, 1)
 
    all_parts_list = []

    # Case 1:  Head + EF2 + Linker2 + EF3 + Tail  (Head and Linker 2 should come from the original LanM sequence where EF2 comes from ; Tail comes from the original LanM sequence where EF3 comes from) (225 cases)
    case1_build = "Head + EF2 + Linker2 + EF3 + Tail".split(" + ")
    all_parts_list += buildParts(selects, case1_build)
    
    # Case 2: Head + EF1 + Linker1 + EF2 + Tail (same logic follows for choosing Head, Linker1, and Tail) (75 cases)
    case2_build = "Head + EF1 + Linker1 + EF2 + Tail".split(" + ")
    all_parts_list += buildParts(selects, case2_build)
    
    # Case 3: Head + EF4 + Linker1 + EF2 + Tail (75 cases)
    case3_build = "Head + EF4 + Linker1 + EF2 + Tail".split(" + ")
    all_parts_list += buildParts(selects, case3_build)
    
    # Case 4: Head + EF1 + Linker1 + EF3 + Tail (75 cases)
    case4_build = "Head + EF1 + Linker1 + EF3 + Tail".split(" + ")
    all_parts_list += buildParts(selects, case4_build)
    
    # Case 5: Head + EF4 + Linker1 + EF3 + Tail (75 cases)
    case5_build = "Head + EF4 + Linker1 + EF3 + Tail".split(" + ")
    all_parts_list += buildParts(selects, case5_build)

    fasta_list = []
    print("Building lanM half design...")
    for i in tqdm(range(len(all_parts_list))):
        fasta_list.append(makeHalfDesign(frags, all_parts_list[i], i))

    writeFastaTupleListToFile(fasta_list, "Half_designs.fasta")
    pass

#_____
# main
def main():
    args = parser()
    combine(args.fasta, args.select)

if __name__ == "__main__":
    main()