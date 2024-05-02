# Author: Yee Chuen Teoh
# python build_ef_data.py
'''
update 5-1-2024:
    - please run this script at the base directory, do not run it in the `python_script` directory.
        i.e. move this file to `/work/ratul1/chuen/Lanmodulin-data`

update 4/16/2024:
    - to build fasta for full lanmodulin and ef hand proteins with pdb name.
'''
import os

def getFileContent(file_name):
    with open(file_name, 'r') as f:
        contents = f.readlines()
        contents = [x.replace("\n", "") for x in contents]
        return contents

def removeElement(given_list, to_remove):
    while to_remove in given_list:
        given_list.remove(to_remove)

proteins = getFileContent("proteins.txt")
ef_locations = getFileContent("EFHands_locations_result_list.txt")
removeElement(proteins, "")
removeElement(ef_locations, "")

#print(proteins[:6])
#print(ef_locations[:6])
ef_locations = [eval(x) for x in ef_locations]

result = []
ef_pdbs = [x for x in os.listdir('ef_pdbs') if '.pdb' in x]
print(len(ef_pdbs))

def getPdbfromList(pdb_id, pdbs):
    file_name = []
    for pdb in pdbs:
        if pdb_id in pdb:
            file_name.append(pdb)
    
    return file_name

mapping = "index_to_tag_batch_no_dict.txt" 
mapping_content = [] 
with open(mapping, 'r') as f:
    contents = f.readlines()
    contents = [x.replace("\n", "") for x in contents if x != "\n"]
    contents = [x.replace("index: ", "") for x in contents]
    contents = [x.split(" (Tag, batch_no) = ") for x in contents]
    contents = [eval(x[1]) for x in contents if x[1] != 'Invalid serial number']
    mapping_content = contents

ef_in_lanM = {}
for i in range(len(ef_locations)):
    idx_p = i * 2
    seq = proteins[idx_p + 1]
    ef_num = 1
    for tuple in ef_locations[i]:
        motif_type = ef_locations[i][tuple][0]

        if len(seq[tuple[0]:tuple[1]]) != 12:
            print(f"EF hand not length of 12: {seq[tuple[0]:tuple[1]]}")
            print(tuple)
            print(ef_locations[i][tuple])
            print(proteins[idx_p])
            print(proteins[idx_p + 1])

        if mapping_content[i][0] not in ef_in_lanM: ef_in_lanM[mapping_content[i][0]] = []
        ef_in_lanM[mapping_content[i][0]].append(f"{tuple}")

        result.append(f"{proteins[idx_p]}_EF{ef_num}|{motif_type}|{tuple}")
        result.append(seq[tuple[0]:tuple[1]])

        ef_num += 1


#print(mapping_content)
#print(proteins)


pdbs = [x for x in os.listdir('pdbs') if '.pdb' in x]
with open('Lanmodulins.fasta', 'w') as f:
    for i in range(0, len(proteins), 2):
        seq = proteins[i + 1]

        mapping = mapping_content[int(i / 2)]
        #print(mapping[0])

        name = getPdbfromList(mapping[0], pdbs)

        if len(name) != 1: raise ValueError("\nFound pdb name not equal to 1")

        ef_portion = "".join(ef_in_lanM[mapping[0]])
        description = f">{proteins[i]}_{mapping[0]}|{ef_portion}|{name[0]}"

        f.write(f"{description}\n")
        f.write(f"{seq}\n")

print(len(mapping_content))
print(len(result) / 8)
if len(mapping_content) != int(len(result) / 8): raise ValueError(f"\nBoth list has different size.\n{len(mapping_content)}, {int(len(result) / 8)}")
#result = [i for i in range(len(result))]

idx_result = 0
for i in range(len(mapping_content)):
    for j in range(4):

        name = getPdbfromList(mapping_content[i][0], ef_pdbs)
        #print(name)
        name = getPdbfromList(f"_EF{j + 1}", name)
        #print(name)
        if len(name) != 1: raise ValueError("\nFound pdb name not equal to 1")

        #print(idx_result)
        result[idx_result] = result[idx_result].replace("_EF", f"_{mapping_content[i][0]}_EF")

        result[idx_result] += f"|{name[0]}"
        #print(result[idx_result])
        idx_result += 2

with open("EFhands.fasta", "w") as f:
    for i, x in enumerate(result):
        if i % 2 == 0:
            f.write(f">{x}\n")
        else:
            f.write(f"{x}\n")

with open("all_proteins.fasta", "w") as f:
    for i in range(0, len(proteins), 2):
        seq = proteins[i + 1]

        mapping = mapping_content[int(i / 2)]
        #print(mapping[0])

        name = getPdbfromList(mapping[0], pdbs)

        if len(name) != 1: raise ValueError("\nFound pdb name not equal to 1")

        ef_portion = "".join(ef_in_lanM[mapping[0]])
        description = f">{proteins[i]}_{mapping[0]}|{ef_portion}|{name[0]}"

        f.write(f"{description}\n")
        f.write(f"{seq}\n")

    for i, x in enumerate(result):
        if i % 2 == 0:
            f.write(f">{x}\n")
        else:
            f.write(f"{x}\n")
