# Author: Yee Chuen Teoh
# python build_ef_struct_data.py
from zipfile import ZipFile
import os
from tqdm import tqdm

# base working directory/path
base_path = "/work/ratul1/chuen/Lanmodulin-data"
 
#Biopython import
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

def extractZip(zip_file):
    with ZipFile(zip_file, 'r') as zip: 
        # printing all the contents of the zip file 
        zip_dir_list = zip.namelist()

        file_to_extract = []
        for file in zip_dir_list:
            if 'rank_005' in file and '.pdb' in file:
                file_to_extract.append(file)
        if len(file_to_extract) != 1:
            print("This zip file has more than 1 rank 005 pdb file.")

        file_extract = file_to_extract[0]
        zip.extract(file_extract)

zip_directory = f"{base_path}/_batches/results" # <-- change path to match your local
batches = os.listdir(zip_directory)
if ".DS_Store" in batches: batches.remove(".DS_Store")
batches.sort(key = lambda x: int(x[-2:]))


pdbs_path = f"{base_path}/pdbs"
if not os.path.exists(pdbs_path):
    os.mkdir(pdbs_path)
os.chdir(pdbs_path)

total = 0
print("Unzipping pdb zip files...")
for batch in tqdm(batches):
    batch_list_dir = os.listdir(f"{zip_directory}/{batch}")
    if ".DS_Store" in batch_list_dir: batch_list_dir.remove(".DS_Store")

    zip_list = []
    for f in batch_list_dir:
        if ".zip" in f:
            zip_list.append(f)
    
    total += len(zip_list)
    #print(len(zip_list), batch)
    for z in zip_list:
        if not os.path.exists(f"{pdbs_path}/{z}"):
            extractZip(f"{zip_directory}/{batch}/{z}") # <-- only extract if it does not exist

print("Total pdbs:", total)


# extract EF hand.
ef_pdb_path = f"{base_path}/pdbs" # <-- change path to match your local
ef_pdb_list = os.listdir(ef_pdb_path)
if ".DS_Store" in ef_pdb_list: ef_pdb_list.remove(".DS_Store")
ef_pdb_list = [f"{ef_pdb_path}/{p}" for p in ef_pdb_list]

def getEFHandPdb(pdb_file, ef_dict, protein, seq):
    ef_set = set()

    pdbparser = PDBParser()
    structure = pdbparser.get_structure('seq', pdb_file)
    chains = {chain.id:seq1(''.join(residue.resname for residue in chain)) for chain in structure.get_chains()}
    query_chain = chains['A']
    
    if seq != query_chain: print("Sequence does NOT match pdb sequence")

    ef_list = list(ef_dict.keys())
    with open(pdb_file, 'r') as f:
        contents = f.readlines()
        #contents = [x.replace("\n", "") for x in contents]

        ef_hands_pdbs = []
        ef = []
        ef_idx = 0
        #print(ef_list[ef_idx][0], ef_list[ef_idx][1])
        for line in contents:
            type = line[:4]
            #chain = line[21:22] # <-- we can ignore chain information since it is not needed
            position = (line[22:26].replace(" ",""))
            if type == "ATOM":
                position = int(position)
                if ef_list[ef_idx][0] + 1 <= position and position < ef_list[ef_idx][1] + 1:
                    #print(position)
                    ef.append(line)

                if position > ef_list[ef_idx][1] + 1:
                    ef_hands_pdbs.append(ef)
                    ef = []
                    ef_idx += 1
                    if ef_idx >= len(ef_list):
                        break

        pdb_name = pdb_file.split("/")[-1][:-4]
        for e in range(len(ef_hands_pdbs)):
            ef = ef_hands_pdbs[e]
            ef_hand_file_name = f"{pdb_name}_EF{e + 1}.pdb"
            with open(ef_hand_file_name, "w") as f:
                f.write(f"REMARK    protein fasta: >id_{protein}\n")
                f.write(f"REMARK    full seq: {seq}\n")
                f.write(f"REMARK    motif type: {ef_dict[ef_list[e]][0]}\n")
                for line in ef:
                    f.write(line)
                f.write(f"END")
            
amino_acid_codes = {
    "ALA": "A","ARG": "R","ASN": "N","ASP": "D","CYS": "C",
    "GLN": "Q","GLU": "E","GLY": "G","HIS": "H","ILE": "I",
    "LEU": "L","LYS": "K","MET": "M","PHE": "F","PRO": "P",
    "SER": "S","THR": "T","TRP": "W","TYR": "Y","VAL": "V"
}

ef_dir = f"{base_path}/ef_pdbs"
if not os.path.exists(ef_dir):
    os.mkdir(ef_dir)
os.chdir(ef_dir)

#print(len(amino_acid_codes)) 
ef_list = [(11, 23), (35, 47), (60, 72), (84, 96)] # <-- TODO: change this.

mapping_file = f"{base_path}/index_to_tag_batch_no_dict.txt"
mapping = {}
value = set()
with open(mapping_file, "r") as f:
    contents = f.readlines()
    contents = [x.replace("\n", "") for x in contents]
    while "" in contents:
        contents.remove("")
    for line in contents:
        equal_seperator = line.split(" = ")
        if equal_seperator[-1] == "Invalid serial number": continue

        #print(equal_seperator)
        tuple = eval(equal_seperator[-1])
        index = equal_seperator[0].split(" (Tag, batch_no)")[0].split("index: ")[-1]

        if index in value: print("duplicate index.")

        mapping[tuple[0]] = int(index)
        value.add(index)

#print(mapping)

def getFileContent(file_name):
    with open(file_name, 'r') as f:
        contents = f.readlines()
        contents = [x.replace("\n", "") for x in contents]
        return contents
def removeElement(given_list, to_remove):
    while to_remove in given_list:
        given_list.remove(to_remove)
ef_locations_file = f"{base_path}/EFHands_locations_result_list.txt"
ef_locations = getFileContent(ef_locations_file)
removeElement(ef_locations, "")
#print(len(ef_locations)) # <-- development and debug purposes

proteins_file = f"{base_path}/proteins.txt"
proteins = getFileContent(proteins_file)
removeElement(proteins, "")
#print(len(proteins)) # <-- development and debug purposes

def getCodeName(pdb_name):
    pdb_name_split = pdb_name.split("_")
    sol = []
    for s in pdb_name_split:
        if s == "rank" or s == "unrelaxed": break
        sol.append(s)

    return "_".join(sol)

print("Scrapping EF pdb from protein...")
total_ef = 0
for i in tqdm(range(len(ef_pdb_list))):
    #if i == 5: break

    pdb_name = ef_pdb_list[i].split("/")[-1][:-4]
    code_key = getCodeName(pdb_name)
    if code_key not in mapping: print(f"Error name {code_key} not in mapping dict")

    idx_map = int(mapping[code_key])

    protein = proteins[idx_map * 2]
    seq = proteins[(idx_map * 2) + 1]

    total_ef += len(list(eval(ef_locations[idx_map]).keys()))
    getEFHandPdb(ef_pdb_list[i], eval(ef_locations[idx_map]), protein, seq)

print("Total ef hands:", total_ef)
