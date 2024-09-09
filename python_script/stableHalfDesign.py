import pandas as pd

excel_file = "/work/ratul1/curwen/Lanmodulin-data/python_script/energy_results1.xlsx"
output_file = "top_stable_energy_results.xlsx"

try:
    df = pd.read_excel(excel_file)
except FileNotFoundError:
    print(f"Error: Excel file not found at {excel_file}")
    exit()

cases = [
    ('ef2', 'ef3'),
    ('ef1', 'ef2'),
    ('ef4', 'ef2'),
    ('ef1', 'ef3'),
    ('ef4', 'ef3')
]

# Function to filter and get top designs for each case
def get_top_designs(df, case, top_n=10):
    case_pattern = rf"\({case[0]},{case[1]}\)"
    filtered_df = df[df['pdb_name'].str.contains(case_pattern, regex=True)]
    top_designs = filtered_df.nsmallest(top_n, 'energy')
    return top_designs

# Get top designs for each case
top_designs_list = []
for case in cases:
    top_designs = get_top_designs(df, case)
    if not top_designs.empty:
        top_designs_list.append(top_designs)

# Combine
if top_designs_list:
    stable_energy_results = pd.concat(top_designs_list)
    print("Combined top designs DataFrame:")
    print(stable_energy_results)

    stable_energy_results.to_excel(output_file, index=False)

    print(f"File containing top entries generated successfully: {output_file}")
else:
    print("No top designs foimport pandas as pd
import matplotlib.pyplot as plt

# File paths
excel_file = "/Users/curwenpeihongtan/Desktop/Work/scatterplot2/top_stable_energy_results.xlsx"
output_file = "top_stable_energy_results.xlsx"

# Load the Excel file into a DataFrame
try:
    df = pd.read_excel(excel_file)
except FileNotFoundError:
    print(f"Error: Excel file not found at {excel_file}")
    exit()

# Define the cases
cases = [
    ('ef2', 'ef3'),
    ('ef1', 'ef2'),
    ('ef4', 'ef2'),
    ('ef1', 'ef3'),
    ('ef4', 'ef3')
]

# Function to filter and get top designs for each case
def get_top_designs(df, case, top_n=10):
    case_pattern = rf"\({case[0]},{case[1]}\)"
    filtered_df = df[df['pdb_name'].str.contains(case_pattern, regex=True)]
    top_designs = filtered_df.nsmallest(top_n, 'energy')
    return top_designs

# Get top designs for each case
top_designs_list = []
for case in cases:
    top_designs = get_top_designs(df, case)
    if not top_designs.empty:
        top_designs_list.append(top_designs)

# Combine all top designs
if top_designs_list:
    stable_energy_results = pd.concat(top_designs_list)
    print("Combined top designs DataFrame:")
    print(stable_energy_results)

    # Save the top designs to a new Excel file
    stable_energy_results.to_excel(output_file, index=False)

    # Print success message
    print(f"File containing top entries generated successfully: {output_file}")

    # Plotting the scatterplot
    plt.figure(figsize=(14, 8))  # Increased figure size for better layout
    for case in cases:
        case_pattern = rf"\({case[0]},{case[1]}\)"
        case_df = stable_energy_results[stable_energy_results['pdb_name'].str.contains(case_pattern, regex=True)]
        plt.scatter(case_df['pdb_name'], case_df['energy'], label=f"{case[0]},{case[1]}", s=30)  # Increased marker size

    plt.xlabel('PDB Name', fontsize=20)  # Increased font size
    plt.ylabel('Energy', fontsize=20)  # Increased font size
    plt.title('Scatterplot of Energy for Different Cases', fontsize=25)  # Increased font size
    plt.xticks(rotation=90, fontsize=10)  # Rotate x-axis labels and increase font size
    plt.yticks(fontsize=15)  # Increase y-axis labels font size
    plt.legend(fontsize=15)  # Increase legend font size
    plt.tight_layout()
    plt.show()

else:
    print("No top designs found for any case.")
und for any case.")
