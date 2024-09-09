import pandas as pd

def load_data(file_path):
    return pd.read_excel(file_path)

# Extract the top 50 stable combinations
def get_top_combinations(df):
    top_combinations = df.nsmallest(50, 'energy')  
    return top_combinations

def save_for_pyrosetta(combinations, output_file):
    combinations.to_csv(output_file, index=False)

def main():
    df = load_data('/work/ratul1/curwen/Lanmodulin-data/python_script/top_stable_energy_results.xlsx')
    top_combinations = get_top_combinations(df)
    save_for_pyrosetta(top_combinations, 'top_combinations.csv')

if __name__ == '__main__':
    main()
