import itertools

# List of top EF hand entries
top_ef1 = ['EF1_seq1', 'EF1_seq2', 'EF1_seq3', 'EF1_seq4', 'EF1_seq5']
top_ef2 = ['EF2_seq1', 'EF2_seq2', 'EF2_seq3', 'EF2_seq4', 'EF2_seq5', 
           'EF2_seq6', 'EF2_seq7', 'EF2_seq8', 'EF2_seq9', 'EF2_seq10', 
           'EF2_seq11', 'EF2_seq12', 'EF2_seq13', 'EF2_seq14', 'EF2_seq15']
top_ef3 = ['EF3_seq1', 'EF3_seq2', 'EF3_seq3', 'EF3_seq4', 'EF3_seq5', 
           'EF3_seq6', 'EF3_seq7', 'EF3_seq8', 'EF3_seq9', 'EF3_seq10', 
           'EF3_seq11', 'EF3_seq12', 'EF3_seq13', 'EF3_seq14', 'EF3_seq15']
top_ef4 = ['EF4_seq1', 'EF4_seq2', 'EF4_seq3', 'EF4_seq4', 'EF4_seq5']

def generate_combinations():
    cases = {}
    
    # Case 1: Head + EF2 + Linker2 + EF3 + Tail
    cases['case_1'] = list(itertools.product(top_ef2, top_ef3))
    
    # Case 2: Head + EF1 + Linker1 + EF2 + Tail
    cases['case_2'] = list(itertools.product(top_ef1, top_ef2))
    
    # Case 3: Head + EF4 + Linker1 + EF2 + Tail
    cases['case_3'] = list(itertools.product(top_ef4, top_ef2))
    
    # Case 4: Head + EF1 + Linker1 + EF3 + Tail
    cases['case_4'] = list(itertools.product(top_ef1, top_ef3))
    
    # Case 5: Head + EF4 + Linker1 + EF3 + Tail
    cases['case_5'] = list(itertools.product(top_ef4, top_ef3))

    return cases

cases = generate_combinations()

for case_name, combinations in cases.items():
    print(f"Combinations for {case_name}:")
    for combo in combinations:
        print(f"Head + {combo[0]} + Linker + {combo[1]} + Tail")
    print("\n")



