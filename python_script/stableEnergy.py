# Author: Curwen Tan Pei Hong

import pandas as pd

excel_file = "ef_hands_updated.xlsx"
df = pd.read_excel(excel_file)

df.columns = ['pdb_name', 'sequences', 'post_relax_energy']

df['pdb_name'] = df['pdb_name'].str.strip()
df['sequences'] = df['sequences'].str.strip()

# Convert 'post_relax_energy' column to numeric,
# Coercion in Python Pandas is the process of converting data from one type to another. 
# This can be done using the to_numeric() function. The to_numeric() function has an errors argument that can be set to coerce
df['post_relax_energy'] = pd.to_numeric(df['post_relax_energy'], errors='coerce')

# Filter the data 
top_ef2 = df[df['pdb_name'].str.contains('EF2.pdb')].nsmallest(15, 'post_relax_energy')

top_ef3 = df[df['pdb_name'].str.contains('EF3.pdb')].nsmallest(15, 'post_relax_energy')

top_ef1 = df[df['pdb_name'].str.contains('EF1.pdb')].nsmallest(5, 'post_relax_energy')

top_ef4 = df[df['pdb_name'].str.contains('EF4.pdb')].nsmallest(5, 'post_relax_energy')

# Combine all filtered data (concatenation)
stableEnergy_result = pd.concat([top_ef2, top_ef3, top_ef1, top_ef4])

output_file = "stableEnergy_result"
stableEnergy_result['pdb_name'].to_csv(output_file, index=False, header=False)

print("File containing top hands generated successfully!")
