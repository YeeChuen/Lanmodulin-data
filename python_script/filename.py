import os

directory = '/work/ratul1/curwen/Lanmodulin-data/ESMFold_results'

# Read original filenames 
with open(os.path.join(directory, 'ESMFold_results_modified'), 'r') as original_filenames_file:
    original_filenames = original_filenames_file.read().splitlines()

# Loop through original filenames and rename the files back to their original names
for original_filename in original_filenames:
    os.rename(os.path.join(directory, original_filename), os.path.join(directory, original_filename))

# Remove "ESMFold_results_modified" file
os.remove(os.path.join(directory, 'ESMFold_results_modified'))


