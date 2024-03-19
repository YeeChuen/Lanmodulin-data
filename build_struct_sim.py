# Author: Yee Chuen Teoh
# python build_struct_sim.py

#______
# imports


#______
# functions
# check for unique seq
def checkUnique():
    with open('EFhands.fasta', 'r') as f:
        contents = f.readlines()
        seq = set()
        for line in contents:
            if line[0] == ">":
                continue
            line = line.replace("\n", "")
            seq.add(line)
        
        print(len(seq))

#______
# main
def main():
    checkUnique()
    pass

if __name__ == "__main__":
    main()
