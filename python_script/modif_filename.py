'''
import os

directory = '/work/ratul1/curwen/Lanmodulin-data/ESMFold_results'

# Create a new file to store modified filenames
with open(os.path.join(directory, 'ESMFold_results_modified'), 'w') as output_file:
    for filename in os.listdir(directory):
        if filename.endswith('.pdb'):
            # Extract the part before the first "|"
            new_filename = filename.split('|')[0].strip() + '.pdb'
            # Rename the file
            os.rename(os.path.join(directory, filename), os.path.join(directory, new_filename))
            # Write the new filename to the new file
            output_file.write(new_filename + '\n')
            '''
'''
import os

# Directory containing the .pdb files
directory = '/work/ratul1/curwen/Lanmodulin-data/ESMFold_results'

# Loop through all .pdb files in the directory
for filename in os.listdir(directory):
    if filename.endswith('.pdb'):
        # Remove extra '.pdb' extensions from the filename
        new_filename = filename.replace('.pdb', '')
        # Rename the file
        os.rename(os.path.join(directory, filename), os.path.join(directory, new_filename))
'''
import os

# Directory containing the .pdb files
directory = '/work/ratul1/curwen/Lanmodulin-data/ESMFold_results'

# Loop through all .pdb files in the directory
for filename in os.listdir(directory):
    if filename.endswith('.pdb'):
        # Remove extra '.pdb' extensions from the filename
        new_filename = filename.replace('.pdb', '')
        # Rename the file
        os.rename(os.path.join(directory, filename), os.path.join(directory, new_filename))
