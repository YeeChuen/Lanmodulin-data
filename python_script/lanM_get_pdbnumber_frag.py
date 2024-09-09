# Author: Yee Chuen Teoh
'''
takes in a pdb file, and fasta file (with build by build_ef_data.py)
where fasta description should contain format like this 
>0_8FNS|(35, 46)(59, 70)(84, 95)(108, 119)|8FNS.pdb
VDIAAFDPDKDGTIDLKEALAAGSAAFDKLDPDKDGTLDAKELKGRVSEADLKKLDPDNDGTLDKKEYLAAVEAQFKAANPDNDGTIDARELASPAGSALVNLIR

as the last two element
index -2: EF portion (inclusive)
index -1: corresponding pdb file

need to return
head, ef, linkers, tail: inclusive range, and length

also assuming pdb only have 1 chain

usage:
    python lanM_get_pdbnumber_frag.py --pdb ../pdbs/8FNS.pdb --fasta ../fasta/wildtype.fasta

update:
'''
#_____
# imports
import argparse
import os

#_____
# functions
def parser():
    parser = argparse.ArgumentParser(description="Read a PDF file from the RAG model.")
    parser.add_argument("--pdb", help=":the protein pdb file", required = True)
    parser.add_argument("--fasta", help=":The fasta file containing the given protein pdb, structure formated from the python script `build_ef_data.py`", required = True)
    return parser.parse_args()

def readFasta(fasta_file):
    with open (fasta_file, 'r') as f:
        contents = f.readlines()
        contents = [x.replace("\n","") for x in contents if x != "\n"]
        return contents

def getEFRange(fasta_context, pdb_path):
    for line in fasta_context:
        if pdb_path.split("/")[-1] == line.split("|")[-1]:
            ef_range_str = line.split("|")[-2]
            ef_range = [eval(f"({r}") for r in ef_range_str.split("(")[1:]]
            return ef_range
    return None

def readPdbNum(pdb_path):
    with open (pdb_path, 'r') as f:
        contents = f.readlines()
        contents = [x.replace("\n","") for x in contents if x != "\n"]
    
    pdb_nums = []
    for line in contents:
        number = int(line[23:32].strip())
        if number not in pdb_nums: pdb_nums.append(number)

    return pdb_nums

def getFragNumbering(pdb_nums, ef_range):
    start = 0

    store = []

    while ef_range:
        ef_tuple = ef_range.pop(0)

        ef_i = pdb_nums.index(ef_tuple[0])
        ef_j = pdb_nums.index(ef_tuple[1])
        print(ef_i, ef_j)
        # 0-2-4-6----11
        add = -17
        store.append([(pdb_nums[start], pdb_nums[ef_i - 1]), (start, ef_i - 1), ef_i - start])
        store.append([(pdb_nums[ef_i], pdb_nums[ef_j]), (ef_i, ef_j ), 
                      [f"{ef_i+0+add}-{ef_i+0+add}", 
                       f"{ef_i+2+add}-{ef_i+2+add}",
                       f"{ef_i+4+add}-{ef_i+4+add}",
                       f"{ef_i+6+add}-{ef_i+6+add}",
                       f"{ef_i+11+add}-{ef_i+11+add}",
                       ],ef_j - ef_i + 1])

        start = ef_j + 1

    store.append([(pdb_nums[start], pdb_nums[len(pdb_nums) - 1]), (start, len(pdb_nums) - 1), len(pdb_nums) - 1 - start + 1])

    total = 0
    for s in store:
        total += s[-1]
        
    return [store, total]


#_____
# main
def main():
    args = parser()
    fasta_context = readFasta(args.fasta)
    # print(fasta_context) # <-- development and debug

    ef_range = getEFRange(fasta_context, args.pdb)
    # print(ef_range) # <-- development and debug

    pdb_nums = readPdbNum(args.pdb)

    [store, total] = getFragNumbering(pdb_nums, ef_range)

    print(len(pdb_nums), total)
    for s in store:
        print(s)
        try:
            print(",".join(s[-2]))
        except:
            continue

    pass

if __name__ == "__main__":
    main()
