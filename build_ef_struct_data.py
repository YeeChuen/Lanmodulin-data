# Author: Yee Chuen Teoh
# python build_ef_struct_data.py
from zipfile import ZipFile
import os

def extractZip(example_zip_file):
    with ZipFile(example_zip_file, 'r') as zip: 
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

example_zip_file = "E:/Files/Coding Software/ChowdhuryVSC/Lanmodulin-data/_batches/results/batch00/ACA18767-1.result.zip"

zip_directory = "E:/Files/Coding Software/ChowdhuryVSC/Lanmodulin-data/_batches/results"
batches = os.listdir(zip_directory)
batches.sort(key = lambda x: int(x[-2:]))

if not os.path.exists("E:/Files/Coding Software/ChowdhuryVSC/Lanmodulin-data/pdbs"):
    os.mkdir("E:/Files/Coding Software/ChowdhuryVSC/Lanmodulin-data/pdbs")
os.chdir("E:/Files/Coding Software/ChowdhuryVSC/Lanmodulin-data/pdbs")

total = 0
for batch in batches:
    batch_list_dir = os.listdir(f"{zip_directory}/{batch}")

    zip_list = []
    for f in batch_list_dir:
        if ".zip" in f:
            zip_list.append(f)
    
    total += len(zip_list)
    print(len(zip_list), batch)
    for z in zip_list:
        if not os.path.exists(f"{zip_directory}/{batch}/{z}"):
            extractZip(f"{zip_directory}/{batch}/{z}") # <-- only extract if it does not exist

# extract EF hand.
pdb_path = "E:/Files/Coding Software/ChowdhuryVSC/Lanmodulin-data/pdbs"
pdb_list = os.listdir(pdb_path)
pdb_list = [f"{pdb_path}/{p}" for p in pdb_list]

def getEFHandPdb(pdb_file, ef_list):
    with open(pdb_file, 'r') as f:
        contents = f.readlines()
        #contents = [x.replace("\n", "") for x in contents]

        ef_hands_pdbs = []
        ef = []
        ef_idx = 0
        print(ef_list[ef_idx][0], ef_list[ef_idx][1])
        for line in contents:
            type = line[:4]
            chain = line[21:22] # <-- we can ignore chain information since it is not needed
            position = (line[22:26].replace(" ",""))
            if type == "ATOM":
                position = int(position)
                if ef_list[ef_idx][0] + 1 <= position and position <= ef_list[ef_idx][1] + 1:
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
            with open(f"{pdb_name}_EF{e + 1}.pdb", "w") as f:
                for line in ef:
                    f.write(line)
            
amino_acid_codes = {
    "ALA": "A","ARG": "R","ASN": "N","ASP": "D","CYS": "C",
    "GLN": "Q","GLU": "E","GLY": "G","HIS": "H","ILE": "I",
    "LEU": "L","LYS": "K","MET": "M","PHE": "F","PRO": "P",
    "SER": "S","THR": "T","TRP": "W","TYR": "Y","VAL": "V"
}

if not os.path.exists("E:/Files/Coding Software/ChowdhuryVSC/Lanmodulin-data/ef_pdbs"):
    os.mkdir("E:/Files/Coding Software/ChowdhuryVSC/Lanmodulin-data/ef_pdbs")
os.chdir("E:/Files/Coding Software/ChowdhuryVSC/Lanmodulin-data/ef_pdbs")

print(len(amino_acid_codes))
ef_list = [(11, 23), (35, 47), (60, 72), (84, 96)] # <-- TODO: change this.

for i in range(len(pdb_list)):
    getEFHandPdb(pdb_list[i], ef_list)
