from pyrosetta import init, pose_from_sequence, get_fa_scorefxn
import pandas as pd

# Initialize PyRosetta
init()

# Load the data prepared for PyRosetta
def load_combinations(file_path):
    return pd.read_csv(file_path)

# Process each combination using PyRosetta
def process_combinations(combinations):
    scorefxn = get_fa_scorefxn()  # Default scoring function
    results = []
    for index, row in combinations.iterrows():
        sequence = row['Combination']  # Ensure the column name matches
        pose = pose_from_sequence(sequence)
        energy = scorefxn(pose)
        results.append((sequence, energy))
    return results

# Main function to control workflow
def main():
    combinations = load_combinations('top_combinations.csv')
    results = process_combinations(combinations)
    for result in results:
        print(f"Combination: {result[0]}, Energy: {result[1]}")

if __name__ == '__main__':
    main()
