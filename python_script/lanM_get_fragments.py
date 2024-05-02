# Author: Yee Chuen Teoh
'''
usage:
    python lanM_get_fragments.py --fasta "/work/ratul1/chuen/Lanmodulin-data/Lanmodulins.fasta"

update:
creation of script (5/1/2024)
    - this script is created to get the 9 fragment of all lanM, as defined:
        A Full LanM protein sequence is arranged into these nine fragments
        Head + EF1 + Linker1 + EF2 + Linker2 + EF3 + Linker3 + EF4 + Tail 
'''
#_____
# imports
import argparse
import os

#_____
# functions
def parser():
    parser = argparse.ArgumentParser(description="Read a PDF file from the RAG model.")
    parser.add_argument("--fasta", help=":The fasta file containing lanM, structure formated from the python script `build_ef_data.py`", required = True)
    return parser.parse_args()

def readFasta(fasta_file):
    with open (fasta_file, 'r') as f:
        contents = f.readlines()
        contents = [x.replace("\n","") for x in contents if x != "\n"]
        return contents

def getFragments(fasta_file):
    fasta_content = readFasta(fasta_file)
    
    fragments = {}
    for i in range(0, len(fasta_content), 2):
        description = fasta_content[i][1:]
        seq = fasta_content[i + 1]

        d_temp = description.split("|")
        lanM_id = "_".join(d_temp[0].split("_")[1:])
        fragments[lanM_id] = {}
        
        fragments[lanM_id]["efs_portion"] = d_temp[1]
        fragments[lanM_id]["file_name"] = d_temp[2]
        
        ep_temp = [f"{x})" for x in d_temp[1].split(")") if x != ""]

        
        non_ef_frags = ['head', 'linker1','linker2','linker3','tail']
        prev = 0
        for j in range(len(ep_temp)):
            tuple = eval(ep_temp[j])
            ef_seq = seq[tuple[0]:tuple[1]]
            fragments[lanM_id][f"ef{j + 1}"] = ef_seq

            non_ef_seq = seq[prev:tuple[0]]
            fragments[lanM_id][non_ef_frags[j]] = non_ef_seq
            prev = tuple[1]

        non_ef_seq = seq[prev:len(seq)]
        fragments[lanM_id][non_ef_frags[j + 1]] = non_ef_seq

        combined_seq = fragments[lanM_id]['head'] + fragments[lanM_id]['ef1'] + fragments[lanM_id]['linker1'] + fragments[lanM_id]['ef2'] + fragments[lanM_id]['linker2'] + fragments[lanM_id]['ef3'] + fragments[lanM_id]['linker3'] + fragments[lanM_id]['ef4'] + fragments[lanM_id]['tail']
        if (seq != combined_seq):
            raise ValueError(f"\nIncorrect fragment split for\n{lanM_id}\nsequence: {seq}\nfragment: {combined_seq}")
        
    return fragments

#_____
# main
def main():
    fragments = getFragments(parser().fasta)
    print(fragments)
    pass

if __name__ == "__main__":
    main()