# Author: Yee Chuen Teoh
# python build_ef_data.py

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

print(proteins[:6])
print(ef_locations[:6])
ef_locations = [eval(x) for x in ef_locations]

result = []
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

        result.append(f"{proteins[idx_p]}_EF{ef_num}|{motif_type}")
        result.append(seq[tuple[0]:tuple[1]])

        ef_num += 1

mapping = "index_to_tag_batch_no_dict.txt" 
mapping_content = [] 
with open(mapping, 'r') as f:
    contents = f.readlines()
    contents = [x.replace("\n", "") for x in contents if x != "\n"]
    contents = [x.replace("index: ", "") for x in contents]
    contents = [x.split(" (Tag, batch_no) = ") for x in contents]
    contents = [eval(x[1]) for x in contents if x[1] != 'Invalid serial number']
    mapping_content = contents

print(len(mapping_content))
print(len(result) / 8)
if len(mapping_content) != int(len(result) / 8): raise ValueError(f"\nBoth list has different size.\n{len(mapping_content)}, {int(len(result) / 8)}")
#result = [i for i in range(len(result))]

idx_result = 0
for i in range(len(mapping_content)):
    for j in range(4):
        #print(idx_result)
        result[idx_result] = result[idx_result].replace("_EF", f"_{mapping_content[i][0]}_EF")
        #print(result[idx_result])
        idx_result += 2

with open("EFhands.fasta", "w") as f:
    for i, x in enumerate(result):
        if i % 2 == 0:
            f.write(f">{x}\n")
        else:
            f.write(f"{x}\n")
