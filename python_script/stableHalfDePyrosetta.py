import os
from pyrosetta import pose_from_pdb, get_fa_scorefxn, init
from pyrosetta.rosetta.core.import_pose import pose_from_pdbstring

init()

pdb_list_file = "pdb_list.txt"
energy_output_file = "energy_results.csv"

def calculate_energy(pdb_file_path):
    try:
        pose = pose_from_pdb(pdb_file_path)
        scorefxn = get_fa_scorefxn()
        energy = scorefxn(pose)
        return energy
    except Exception as e:
        print(f"Error processing file {pdb_file_path}: {e}")
        return None

if os.path.exists(pdb_list_file):
    with open(pdb_list_file, 'r') as f:
        pdb_files = [line.strip() for line in f.readlines()]
    
    energies = []
    for pdb_file in pdb_files:
        energy = calculate_energy(pdb_file)
        if energy is not None:
            energies.append((os.path.basename(pdb_file), energy))
    
    # Save the energies to a CSV file
    with open(energy_output_file, 'w') as f:
        f.write("pdb_name,energy\n")
        for pdb_name, energy in energies:
            f.write(f"{pdb_name},{energy}\n")
    
    print(f"Energy results saved to {energy_output_file}")
else:
    print(f"PDB list file not found: {pdb_list_file}")
