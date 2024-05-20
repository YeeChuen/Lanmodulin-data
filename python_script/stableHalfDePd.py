import pandas as pd
import os

fasta_file = "/work/ratul1/curwen/Lanmodulin-data/fasta/Half_designs.fasta"
pdb_dir = "/work/ratul1/curwen/Lanmodulin-data/ESMFold_results"
pdb_list_file = "pdb_list.txt"

# Check if the FASTA file exists
if not os.path.exists(fasta_file):
    print(f"File not found: {fasta_file}")
else:
    print(f"File found: {fasta_file}")

def parse_fasta(file_path):
    records = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line[1:].strip()
                parts = header.split('|')
                pdb_name = parts[0]
                sequences = next(f).strip()
                records.append((pdb_name, sequences))
    return records

# Read the FASTA file and parse records
if os.path.exists(fasta_file):
    records = parse_fasta(fasta_file)
    print(f"Total records parsed: {len(records)}")

    # Create a DataFrame
    df = pd.DataFrame(records, columns=['pdb_name', 'sequences'])
    print("DataFrame created:")
    print(df.head())

    df['pdb_name'] = df['pdb_name'].str.strip()
    df['sequences'] = df['sequences'].str.strip()

    pdb_list = []
    for index, row in df.iterrows():
        pdb_file_path = os.path.join(pdb_dir, f"{row['pdb_name']}.pdb")
        if os.path.exists(pdb_file_path):
            pdb_list.append(pdb_file_path)
        else:
            print(f"PDB file not found: {pdb_file_path}")

    # Save to a file
    with open(pdb_list_file, 'w') as f:
        for pdb_file in pdb_list:
            f.write(pdb_file + "\n")

    print(f"PDB list saved to {pdb_list_file}")
else:
    print("FASTA file not found. Please check the file path.")
