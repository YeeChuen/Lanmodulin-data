import pandas as pd
from Bio.PDB import PDBParser

def parse_pdb_file(file_path):
    parser = PDBParser()
    structure = parser.get_structure('PDB', file_path)
    atom_details = []

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom_details.append({
                        "atom_number": atom.serial_number,
                        "atom_name": atom.name,
                        "residue_name": residue.resname,
                        "chain_id": chain.id,
                        "residue_number": residue.id[1],
                        "x": atom.coord[0],
                        "y": atom.coord[1],
                        "z": atom.coord[2],
                        "occupancy": atom.occupancy,
                        "temp_factor": atom.bfactor,
                    })
    return atom_details

def process_pdb_files_from_list(file_list_path):
    # Process a list of PDB files from a CSV and collect atom details.
    df_files = pd.read_csv(file_list_path)
    all_atom_details = []
    for index, row in df_files.iterrows():
        # Construct the full path to each PDB file, assuming they are in the same directory as the CSV file.
        pdb_file_path = file_list_path.rsplit('/', 1)[0] + '/' + row['pdb_name']
        atom_details = parse_pdb_file(pdb_file_path)
        all_atom_details.extend(atom_details)
    return all_atom_details

def main():
    file_list_path = 'ratul1/curwen/Lanmodulin-data/python_script/RF_diffusion job/top_combinations.csv'

    results = process_pdb_files_from_list(file_list_path)
    
    # Create a DataFrame from the results and save to CSV
    df_results = pd.DataFrame(results)
    df_results.to_csv('/work/ratul1/curwen/Lanmodulin-data/python_script/RF_diffusion job/extracted_atom_details.csv', index=False)
    print(df_results.head())  

if __name__ == '__main__':
    main()
