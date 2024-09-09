import itertools
from pyrosetta import init, pose_from_sequence, get_fa_scorefxn

# Initialize PyRosetta
init()

# Define parts for combinations
heads = ['Head']
ef2_seqs = ['EF2_seq1', 'EF2_seq2', 'EF2_seq3', 'EF2_seq4', 'EF2_seq5']
ef3_seqs = ['EF3_seq1', 'EF3_seq2', 'EF3_seq3', 'EF3_seq4', 'EF3_seq5']
linkers = ['Linker']
tails = ['Tail']

combinations = list(itertools.product(heads, ef2_seqs, linkers, ef3_seqs, tails))

# Write combinations to a file
with open('protein_combinations.txt', 'w') as file:
    for combo in combinations:
        sequence = ' + '.join(combo) + '\n'
        file.write(sequence)

print(f"Total combinations written: {len(combinations)}")

def process_combinations(filename):
    scorefxn = get_fa_scorefxn()  # Full-atom score function for energy calculations
    with open(filename, 'r') as file:
        for line in file:
            sequence = line.strip().replace(' + ', '')
            pose = pose_from_sequence(sequence)
            energy = scorefxn(pose)
            print(f"Sequence: {sequence}, Energy: {energy}")

