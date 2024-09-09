import os
import csv
import subprocess

# Define paths
sequences_path = "/work/ratul1/curwen/Lanmodulin-data/ESMFold_results"
csv_file_path = "top_combinations.csv"
output_base_path = "/work/ratul1/curwen/Lanmodulin-data/python_script/output_file"

# Function to read and parse the CSV file
def parse_csv(csv_file_path):
    pdb_entries = []
    with open(csv_file_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            pdb_entries.append((row['pdb_name'], float(row['energy'])))
    return pdb_entries

# Function to read and parse the sequences
def parse_sequences(sequences_path):
    sequences = {}
    for filename in os.listdir(sequences_path):
        if filename.endswith(".pdb"):
            with open(os.path.join(sequences_path, filename)) as file:
                sequences[filename] = file.read()
    return sequences

# Function to extract EF-hand ranges from the PDB filename
def extract_ef_ranges(pdb_name):
    ef_part = pdb_name.split('_')[1]
    ef_ranges = ef_part.replace('ef', '').replace('(', '').replace(')', '').replace(',', ' ')
    return ef_ranges

# Function to run RFDiffusion
def run_rfdiffusion(pdb_file_path, output_dir, sequence_name, ef_ranges):
    # Construct the RFDiffusion command with constraints to keep EF-hand sequences fixed
    rfdiffusion_command = [
        "rfdiffusion",
        f"--input_pdb={pdb_file_path}",
        f"--output_dir={output_dir}",
        f"--sequence_name={sequence_name}",
        f"--ef_ranges={ef_ranges}"
    ]
    
    # Using subprocess.run with a list to avoid shell injection issues
    subprocess.run(rfdiffusion_command)

# Main function to run the entire process
def main():
    pdb_entries = parse_csv(csv_file_path)
    sequences = parse_sequences(sequences_path)

    for i, (pdb_name, energy) in enumerate(pdb_entries):
        output_dir = os.path.join(output_base_path, str(i))
        os.makedirs(output_dir, exist_ok=True)
        
        pdb_content = sequences.get(pdb_name)
        if pdb_content:
            pdb_file_path = os.path.join(output_dir, pdb_name)
            with open(pdb_file_path, 'w') as pdb_file:
                pdb_file.write(pdb_content)
            
            ef_ranges = extract_ef_ranges(pdb_name)
            run_rfdiffusion(pdb_file_path, output_dir, pdb_name, ef_ranges)
            # Add logic here to process RFDiffusion output and ensure EF-hand sequences remain fixed

if __name__ == "__main__":
    main()
