# Author: Yee Chuen Teoh
'''
Usage: 
    python report_ef_data.py


# use SE3nv.
    pip3 -q install git+https://github.com/sokrypton/ColabDesign.git@v1.1.1
'''

#______
# imports
import json
from tqdm import tqdm
import re
import os

# ProteinMPNN
from colabdesign.mpnn import mk_mpnn_model

#______
# functions
def getDictFromFile(file):
    sol = {}
    with open(file, 'r') as f:
        contents = f.readlines()
        contents = [x.replace("\n", "") for x in contents if x != "\n"]
        # print(len(contents))
        print(f"Loading from {file}...")
        for x in tqdm(contents):
            sol[eval(x.split(":")[0])] = float(x.split(":")[1])
            
    return sol

def readFastaFile(file):
    with open(file, 'r') as f:
        contents = f.readlines()
        contents = [x.replace("\n", "") for x in contents if x != "\n"]
        return contents 

def getIdMapping():
    protein_mapping = {}
    protein_file = "proteins.txt"
    with open(protein_file, 'r') as f:
        contents = f.readlines()
        contents = [x.replace("\n", "") for x in contents if x != "\n"]
        for i in range(0, len(contents), 2):
            if contents[i] in protein_mapping: print("WARNING, protein has been added to map.")
            protein_mapping[contents[i]] = contents[i + 1]
        #print(protein_mapping)
    
    mapping = "index_to_tag_batch_no_dict.txt" 
    mapping_content = [] 
    mapping_dict = {}
    with open(mapping, 'r') as f:
        contents = f.readlines()
        contents = [x.replace("\n", "") for x in contents if x != "\n"]
        contents = [x.replace("index: ", "") for x in contents]
        contents = [x.split(" (Tag, batch_no) = ") for x in contents]
        for x in contents:
            if x[1] != 'Invalid serial number':
                mapping_dict[eval(x[1])[0]] = protein_mapping[x[0]]
        #print(mapping_dict)
    return mapping_dict, protein_mapping

def getMotifInfo():
    AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    X = r'[ARNDCQEGHILKMFPSTWYV]'
    Y = r'[ARNCQEGHILKMFPSTWYV]'
    DTE = r'[DTE]'

    motif_info = {
        "motif_dA": rf'[DN]-*{X}-*[DN]-*{X}-*D-*{X}-*{X}-*{X}-*{Y}-*{X}-*{Y}-*E',
        "motif_dB" : rf'[DN]-*{X}-*D-*{X}-*{Y}-*{X}-*{X}-*{X}-*D-*{X}-*{Y}-*E',
        "motif_dC" : rf'[DN]-*{X}-*{Y}-*{X}-*D-*{X}-*{X}-*{X}-*D-*{X}-*{Y}-*E',
        "motif_dD" : rf'[DN]-*{X}-*{Y}-*{X}-*D-*{X}-*{X}-*{X}-*{Y}-*{X}-*D-*E',
        "motifE" : rf'E-*{X}-*D-*{X}-*D-*{X}-*{X}-*{X}-*{DTE}-*{X}-*{X}-*E',
        "motif_and_AB" : rf'[DN]-*{X}-*D-*{X}-*D-*{X}-*{X}-*{X}-*D-*{X}-*{Y}-*E',
        "motif_and_AD" : rf'[DN]-*{X}-*D-*{X}-*D-*{X}-*{X}-*{X}-*{Y}-*{X}-*D-*E',
        "motif_and_BD" : rf'[DN]-*{X}-*D-*{X}-*D-*{X}-*{X}-*{X}-*D-*{X}-*D-*E',
        "motif_100_248_413" : rf'N-*[LVR]-*I-*[RK]-*[GR]-*K-*G-*L-*D-*R-*[IV]-*E',
        "motif_100" : rf'D-*R-*[R]-*[GR]-*D-*G-*K-*L-*N-*[IV]-*F-*E',
        "motif_215" : rf'G-*K-*K-*P-*T-*D-*T-*L-*T-*V-*Q-*E',
        "motif_428" : rf'D-*L-*D-*K-*R-*G-*K-*L-*W-*L-*G-*E',
        "motif_246" : rf'D-*[AILHK]-*[IR]-*K-*[DR]-*G-*[KV]-*[LI]-*[TSW]-*[AEGLV]-*[AQKVG]-*E',
        "motif_130" : rf'N-*P-*G-*H-*K-*D-*N-*L-*T-*R-*D-*Q',
        "motif_253_565" : rf'Q-*K-*D-*Q-*D-*D-*T-*L-*D-*[PR]-*K-*E',
        "motif_236" : rf'Y-*T-*D-*S-*D-*G-*T-*L-*D-*H-*I-*E',
        "motif_342" : rf'A-*G-*G-*G-*P-*D-*K-*T-*I-*K-*I-*E',
        "motif_353" : rf'D-*L-*I-*K-*G-*R-*G-*I-*S-*L-*G-*E',
        "motif_432" : rf'D-*R-*G-*K-*D-*G-*T-*L-*T-*A-*R-*E',
        "motif_548" : rf'K-*L-*F-*D-*N-*G-*T-*L-*D-*L-*A-*E'
    }

    motif_info_fix = {
        "motif_dA":             rf'_?_?_??????_',
        "motif_dB" :            rf'_?_?????_??_',
        "motif_dC" :            rf'_???_???_??_',
        "motif_dD" :            rf'_???_?????__',
        "motifE" :              rf'_?_?_??????_',
        "motif_and_AB" :        rf'_?_?_???_??_',
        "motif_and_AD" :        rf'_?_?_?????__',
        "motif_and_BD" :        rf'_?_?_???_?__',
        "motif_100_248_413" :   rf'_?_??_____?_',
        "motif_100" :           rf'___?_____?__',
        "motif_215" :           rf'____________',
        "motif_428" :           rf'____________',
        "motif_246" :           rf'_??_?_?????_',
        "motif_130" :           rf'____________',
        "motif_253_565" :       rf'_________?__',
        "motif_236" :           rf'____________',
        "motif_342" :           rf'____________',
        "motif_353" :           rf'____________',
        "motif_432" :           rf'____________',
        "motif_548" :           rf'____________'
    }

    return motif_info, motif_info_fix

    
def reportRMSDPairInfo():
    print("--------------------------")
    ef_rmsd_file = "ef_rmsd_result_v2.txt"
    ef_gr2_rmsd_file = "ef_equal_seq_rmsd_gr2_v2.txt"
    ef_le2_rmsd_file = "ef_equal_seq_rmsd_le2_v2.txt"

    rmsd_dict = getDictFromFile(ef_rmsd_file)
    rmsd_gr2_dict = getDictFromFile(ef_gr2_rmsd_file)
    rmsd_le_2dict = getDictFromFile(ef_le2_rmsd_file)

    print(f"Total ef hand pairs: {len(rmsd_dict)}")
    print(f"Ef hand equal sequence, RMSD > 0.2: {len(rmsd_gr2_dict)}")
    print(f"Ef hand equal sequence, RMSD <= 0.2: {len(rmsd_le_2dict)}")

def reportCystine():
    print("--------------------------")
    ef_hand_file_fasta = "EFhands.fasta"
    fasta_list = readFastaFile(ef_hand_file_fasta)

    ef_with_cys = []
    ef_w_protein_with_cys = []
 
    mapping_dict, protein_mapping = getIdMapping()
    
    pid_set = set()
    for i in range(0, len(fasta_list), 2):
        if "C" in fasta_list[i + 1]:
            #print(f"This sequence has C:\n{fasta_list[i]}\n{fasta_list[i + 1]}")
            ef_with_cys.append(fasta_list[i])
        else:
            pid = "_".join(fasta_list[i][1:].split("|")[0].split("_")[1:-1])
            if pid not in mapping_dict:
                raise ValueError(f"\nProtein id {pid} is not present in mapping content dictionary.\nPotential issue in scrapping the pid.\n{fasta_list[i]}")
            pid_set.add(pid)
            if "C" in mapping_dict[pid]:
                ef_w_protein_with_cys.append(fasta_list[i])


    if len(pid_set) != len(protein_mapping): raise ValueError("\npid_set not equal size as protein_mapping.")

    print(f"Total Number of EF with CYS: {len(ef_with_cys)}")
    print(f"Total Number of EF w/o CYS, but protein contain CYS: {len(ef_w_protein_with_cys)}")
    print(f"Total Number of EF hands: {int(len(fasta_list) / 2)}")

def reportUniqueSeqInMotif():
    print("--------------------------")
    ef_hand_file_fasta = "EFhands.fasta"
    fasta_list = readFastaFile(ef_hand_file_fasta)

    motif_dict = {}
    for i in range(0, len(fasta_list), 2):
        motif_type = fasta_list[i].split("|")[-1]
        if motif_type not in motif_dict: motif_dict[motif_type] = set()
        motif_dict[motif_type].add(fasta_list[i + 1])
    
    print("Number of unique sequence in each motif:")
    for motif in motif_dict:
        print(f"{motif}: {len(motif_dict[motif])}")

# fix = "1-54,60-80" <-- in this format, the position number reflect number in pdb file.
def proteinMPNN(pdb_file, fix, batch = 10, temperature = 0.1):
    mpnn_model = mk_mpnn_model()
    mpnn_model.prep_inputs(pdb_filename=pdb_file, fix_pos=fix, inverse = True)
    # Suggested values 0.1, 0.15, 0.2, 0.25, 0.3.
    samples = mpnn_model.sample_parallel(batch = batch, temperature=temperature)
    sample_seqs = samples['seq']
    return sample_seqs

def getFixPosition(pos = (3, 15), motif_seq = "?_????__????"):
    if (pos[1] - pos[0]) != len(motif_seq): raise ValueError(f"\nMismatch position range and string length,\n{pos}\n{motif_seq}\n{(pos[1] - pos[0])}, {len(motif_seq)}")

    number = []
    for i in range(len(motif_seq)):
        if motif_seq[i] == "?": number.append(pos[0] + i)
    
    return generate_range_string(number)

def generate_range_string(lst):
    ranges = []
    start = lst[0]
    end = lst[0]
    for num in lst[1:]:
        if num == end + 1:
            end = num
        else:
            if start == end:
                ranges.append(f"{start}-{end}")
            else:
                ranges.append(f"{start}-{end}")
            start = num
            end = num
    if start == end:
        ranges.append(f"{start}-{end}")
    else:
        ranges.append(f"{start}-{end}")
    return ','.join(ranges)

def motifMPNN():
    print("--------------------------")
    ef_hand_file_fasta = "EFhands.fasta"
    fasta_list = readFastaFile(ef_hand_file_fasta)

    motif_dict = {}
    for i in range(0, len(fasta_list), 2):
        motif_type = fasta_list[i].split("|")[-1]
        if motif_type not in motif_dict: motif_dict[motif_type] = []
        
        pid = "_".join(fasta_list[i][1:].split("|")[0].split("_")[1:-1]) # <-- proven to work from reportCystine(), else it will throw error during reportCystine()
        motif_dict[motif_type].append((pid, fasta_list[i + 1]))

    mapping_dict, protein_mapping = getIdMapping()
    motif_info, motif_info_fix = getMotifInfo()

    #print(mapping_dict) # {'WP_203154092-1': 'EAKPKTDRTIAAVDTDSDGTIDLAEV'}
    #print(protein_mapping) # {'621': 'QKSAVITAFDPDKD}'

    redesign = []
    motif_le2 = []

    for motif_type in tqdm(list(motif_dict.keys())):
        pattern = re.compile(motif_info[motif_type])

        ef_list_in_motif = motif_dict[motif_type]
        if len(ef_list_in_motif) < 3: 
            motif_le2.append(motif_type)
            continue # only check motif that has 3 or more sequence.
        #print(ef_list_in_motif)
        sequences = set([tuple[1] for tuple in ef_list_in_motif])
        #print(sequences)
        
        prev_max_count = 0
        max_pid = ""

        for tuple in ef_list_in_motif:
            if not pattern.match(tuple[1]): raise ValueError(f"\nEF does not match motif,\n{motif_type}\n{tuple[1]}")

            protein_seq = mapping_dict[tuple[0]]
            #print(protein_seq)
            count = 0

            for seq in sequences:
                if seq in protein_seq:
                    #print(f"{tuple[0]}\n{seq}\n{protein_seq}") # <-- for checker.
                    count += 1
            if count > prev_max_count:
                prev_max_count = count
                max_pid = tuple[0]
        
        '''
        pattern = re.compile(motif_dict["motif_dA"])
        print(pattern.match("DPDDDKTLTKEE"))
        print(pattern.match("EKDNDGTVDRKE"))
        '''
        full_seq = mapping_dict[max_pid]

        number_ef = []
        for i in range(len(full_seq)):
            portion = full_seq[i:i + 12] # <-- 12 because all ef length is 12
            if len(portion) < 12: continue
            if pattern.match(portion): number_ef.append((i + 1,i + 12 + 1))
        
        fix_list = []
        for pos_range in number_ef:
            fix_list.append(getFixPosition(pos = pos_range, motif_seq = motif_info_fix[motif_type]))
        pos_fix = ",".join(fix_list)

        protein_list = os.listdir('pdbs')
        #print(protein_list)
        pdb_file = ""
        for x in protein_list:
            if max_pid in x:
                if pdb_file: raise ValueError(f"\nFound more than 1 matching pdb with similar pdbid,\n{max_pid}")
                pdb_file = f"pdbs/{x}"

        mpnn_seq = proteinMPNN(pdb_file, pos_fix, batch = 10, temperature = 0.2)
        redesign.append(f"\n--- motif: {motif_type} ---")
        redesign.append(f"Max same motif in protein: {max_pid}")
        redesign.append(f"Number of same motif in protein: {prev_max_count}")
        redesign.append(f"Ef hand portion: {number_ef}")
        redesign.append(f"Original protein sequence:\n    {full_seq}")
        redesign.append(f"Fix range used for MPNN:\n    {pos_fix}")
        mpnn_str = "\n    ".join(mpnn_seq)
        redesign.append(f"MPNN designed protein sequence:\n    {mpnn_str}")

    with open("EF_MPNN_design.txt", "w") as f:
        for line in redesign:
            f.write(f"{line}\n")
            print(line)

    print(f"\nMotif with less than 3 element:\n{motif_le2}")
    print(f"No. Motif with less than 3 element:\n{len(motif_le2)}")



#______
# main
def main():
    motifMPNN()
    reportCystine()
    reportUniqueSeqInMotif()
    reportRMSDPairInfo()



if __name__ == "__main__":
    main()
