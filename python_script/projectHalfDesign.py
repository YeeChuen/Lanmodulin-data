import os

fasta_directory = "/work/ratul1/curwen/Lanmodulin-data/fasta" # <-- this is the mistake, this directory is your project directory not fasta directory

output_directory = "/work/ratul1/curwen/Lanmodulin-data/ESMFold_results"

def run_esmFold(fasta_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    esmFold_command = f'esmfold_inference.py -i {fasta_file}  -o {output_dir}' # <-- this part has error but we're not debugging it right now
    print("ESM Fold Command:", esmFold_command) 
    os.system(esmFold_command)

def main():
    half_design_file = os.path.join(fasta_directory, "Half_designs.fasta")
    if os.path.exists(half_design_file):
        print(half_design_file)
        run_esmFold(half_design_file, output_directory)
    else:
        print("Half_designs.fasta not found in the directory.")

if __name__ == "__main__":
    main()
