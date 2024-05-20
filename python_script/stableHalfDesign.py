import pandas as pd

excel_file = "/Users/curwenpeihongtan/Desktop/Work/scatterplot2/energy_results1.xlsx"
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
    print("No top designs found for any case.")
