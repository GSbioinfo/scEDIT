import pandas as pd
import numpy as np
import re
from Bio.SeqUtils import MeltingTemp as mt
from Bio import pairwise2
from itertools import product

def reverse_complement(seq):
    return str(seq.translate(str.maketrans("ACGT", "TGCA"))[::-1])

def one_hot_encode(seq, length=30):
    """One-hot encode a mononucleotide sequence (e.g. 30-mer) as (length, 4)"""
    base_order = ['A', 'C', 'G', 'T']
    encoding = np.zeros((length, 4), dtype=int)
    
    seq = seq.upper()
    if len(seq) != length:
        raise ValueError(f"Expected sequence of length {length}, got {len(seq)}.")

    for i, base in enumerate(seq):
        if base in base_order:
            encoding[i, base_order.index(base)] = 1
        # You could handle N or other IUPAC codes here if needed
    return encoding
def one_hot_encode_dinuc(seq):
    """One-hot encode a 30-mer as (29, 16) dinucleotide matrix"""
    dinucs = [''.join(p) for p in product('ACGT', repeat=2)]
    encoding = np.zeros((len(seq) - 1, 16), dtype=int)

    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2].upper()
        if dinuc in dinucs:
            encoding[i, dinucs.index(dinuc)] = 1
    return encoding


def get_dinucleotide_features(seq):
    """Generate dinucleotide features (positional counts)"""
    seq = seq.upper()
    features = {}
    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        features[f'dinuc_{i}_{dinuc}'] = 1
    return features

def compute_melting_temp(seq):
    try:
        return mt.Tm_NN(seq)
    except:
        return 0.0

def count_base_edits(original, edited):
    count = 0
    for a, b in zip(original.upper(), edited.upper()):
        if a != b and a in 'ACGT' and b in 'ACGT':
            count += 1
    return count

def extract_features(row):
    features = {}

    guide = row['guide'].upper()
    context = row['context'].upper()
    pam = row['pam'].upper()
    original = row['original_target'].upper()
    edited = row['edited_target'].upper()

    # Validate lengths
    if len(guide) != 20 or len(original) != 24 or len(edited) != 24:
        return None

    # One-hot encodings of mono and di
    features.update({
    f'OH_{i}_{b}': v
    for i, row in enumerate(one_hot_encode(context))
    for b, v in zip(['A', 'C', 'G', 'T'], row)
    })

    features.update({
        f'DI_{i}_{d}': v
        for i, row in enumerate(one_hot_encode_dinuc(context))
        for d, v in zip([''.join(p) for p in product('ACGT', repeat=2)], row)
    })

    # GC content
    gc_content = (guide.count('G') + guide.count('C')) / len(guide)
    features['gc_content'] = gc_content

    # Melting temp
    features['Tm_guide'] = compute_melting_temp(guide)
    features['Tm_target'] = compute_melting_temp(original)

    # PAM features
    features['pam_' + pam] = 1

    # Dinucleotide features
    features.update(get_dinucleotide_features(context))

    # Position-specific nucleotides in guide
    for i, base in enumerate(guide):
        features[f'pos{i}_{base}'] = 1

    # Base edit count
    features['num_edits'] = row['mismatches']

    return features

def extract_doench_features(guide,seq_30mer):
    seq_30mer = seq_30mer.upper()
    features = {}

    if len(seq_30mer) != 30:
        raise ValueError("Input sequence must be exactly 30 nt")
    if len(guide) != 20:
        raise ValueError("Input guide sequence must be exactly 20 nt")

    #guide = seq_30mer[2:22]
    # One-hot encodings for mon and di 
    features.update({
    f'OH_{i}_{b}': v
    for i, row in enumerate(one_hot_encode(seq_30mer))
    for b, v in zip(['A', 'C', 'G', 'T'], row)
    })

    features.update({
        f'DI_{i}_{d}': v
        for i, row in enumerate(one_hot_encode_dinuc(seq_30mer))
        for d, v in zip([''.join(p) for p in product('ACGT', repeat=2)], row)
    })
    for i, nt in enumerate(seq_30mer):
        for base in "ACGT":
            features[f"pos_{i}_{base}"] = int(nt == base)

    for i in range(len(guide) - 1):
        dinuc = guide[i:i+2]
        features[f"dinuc_{i}_{dinuc}"] = 1

    gc_count = guide.count('G') + guide.count('C')
    features["gc_count"] = gc_count
    features["gc_low"] = int(gc_count < 10)
    features["gc_high"] = int(gc_count > 10)

    features["TTTT_present"] = int("TTTT" in guide)

    features["Tm_guide"] = mt.Tm_NN(guide)
    features["Tm_5prime"] = mt.Tm_NN(guide[:10])
    features["Tm_3prime"] = mt.Tm_NN(guide[-10:])

    for base in "ACGT":
        features[f"{base}_count"] = guide.count(base)

    features["G_at_20"] = int(guide[19] == "G")
    features["C_at_16"] = int(guide[15] == "C")
    features["A_at_3"] = int(guide[2] == "A")
    features["A_at_4"] = int(guide[3] == "A")
    features["G_at_18"] = int(guide[17] == "G")
    features["GT_19_20"] = int(guide[18:20] == "GT")
    # Base edit count
    
    return features


def compare_seqs(seq1, seq2):
    alignment = pairwise2.align.localms(seq1, seq2, 2, -1, -2, -2)[0]
    aligned1, aligned2 = alignment.seqA, alignment.seqB
    
    mismatches = 0
    indels = 0

    for a, b in zip(aligned1, aligned2):
        if a != b:
            if a == '-' or b == '-':
                indels += 1
            else:
                mismatches += 1
    if(indels>0):
        identity = -(1-sum(a == b for a, b in zip(aligned1, aligned2)) / len(aligned1))
    else:
        identity = (1-sum(a == b for a, b in zip(aligned1, aligned2)) / len(aligned1))
    #identity = sum(a == b for a, b in zip(aligned1, aligned2)) / len(aligned1)
    return identity, mismatches, indels, aligned1, aligned2

def extract_guide(guide,context_seq):
    gstart= context_seq.find(guide)
    gend = gstart+len(guide)
    return context_seq[0:gstart]+"-"+context_seq[gend:]

def process_input(input_csv, output_csv='doench_features.csv'):
    df = pd.read_csv(input_csv)
    all_features = []
    for _, row in df.iterrows():
        feats = extract_features(row)
        if feats:
            feats['guide'] = row['guide']
            feats['edited_target'] = row['edited_target']
            feats['original_target'] = row['original_target']
            all_features.append(feats)

    df_out = pd.DataFrame(all_features).fillna(0)
    df_out.to_csv(output_csv, index=False)
    print(f"[✔] Features saved to {output_csv}")

# --------------------
# Example usage
# --------------------
if __name__ == "__main__":
    #process_input("input_sequences.csv")
    output_csv = "OneHot_encoded_featture_input_GATA1HbF1.csv"
    all_aiml_data = pd.read_csv("/home/gajendra/Dropbox/Tapestri_workflow/PRJNA889818/GATA1_aiml_data_set_from_all_replicate.csv") 
    all_aiml_data = pd.concat([all_aiml_data,pd.read_csv("/home/gajendra/Dropbox/Tapestri_workflow/PRJNA889818/HbF1_aiml_data_set_from_all_replicate.csv")])
    #all_aiml_data = pd.read_csv("/home/gajendra/Dropbox/Tapestri_workflow/PRJNA889818/HbF1_aiml_data_set_from_all_replicate.csv")
    #grna_contex_pam = pd.read_csv("HBF1_aiml_base_editor_targets.csv")
    all_aiml_data[["Edit_Efficacy", "mismatches", "indels", "aligned_orig", "aligned_edit"]] = all_aiml_data.apply(
        lambda row: pd.Series(compare_seqs(row["REFSEQ"], row["SEQ"])), axis=1
    )
    grna_contex_pam = pd.read_csv("GATA1_aiml_base_editor_targets.csv")
    grna_contex_pam = pd.concat([grna_contex_pam,pd.read_csv("HBF1_aiml_base_editor_targets.csv")])
    grna_contex_pam["gRNASeq"] = grna_contex_pam["gRNA"]
    #grna_contex_pam["gRNAinContext"] = grna_contex_pam.apply(lambda row: row["gRNASeq"].upper() not in row["context_30nt"].upper(), axis=1)


    #grna_contex_pam_filtered["gRNAinContext"] = grna_contex_pam_filtered.apply(lambda row: row["gRNASeq"].upper() in row["context_30nt"].upper(), axis=1)
    new_contex_aiml_data = pd.merge(grna_contex_pam, all_aiml_data, on ="gRNASeq", how='outer')
    missed_gRNAS = new_contex_aiml_data.loc[new_contex_aiml_data["gRNA"].isna()].copy(deep=True)

    matched_gRNAS = new_contex_aiml_data.loc[~new_contex_aiml_data["gRNA"].isna()].copy(deep=True)

    #matched_gRNAS[["Edit_Efficacy", "mismatches", "indels", "aligned_orig", "aligned_edit"]] = matched_gRNAS.apply(
    #    lambda row: pd.Series(compare_seqs(row["REFSEQ"], row["SEQ"])), axis=1
    #)

    edited_df = matched_gRNAS.loc[matched_gRNAS['EDTISTATUS']=='Y'].copy(deep=True)


    unedited_df = matched_gRNAS.loc[matched_gRNAS['EDTISTATUS']=='N'].copy(deep=True)

    matched_gRNAS = matched_gRNAS.reset_index(drop=True)
    #matched_gRNAS['Newguide'] = matched_gRNAS[['gRNA','context_30nt']].apply(lambda row: extract_guide(row['gRNA'],row['context_30nt']), axis=1)
    matched_gRNAS['guide_id'] = matched_gRNAS['gRNA']+["_id_" + str(i) for i in matched_gRNAS.index]
    all_features = []
    matched_gRNAS_out =pd.DataFrame()
    for _, row in matched_gRNAS.iterrows():
        feats = extract_doench_features(row['gRNA'],row['context_30nt'])
        if feats:
            feats['guide_id'] = row['guide_id']
            feats['guide'] = row['gRNA']
            #feats['edited_target'] = row['SEQ']
            #feats['original_target'] = row['REFSEQ']
            feats['editsStatus']=row['EDTISTATUS']
            feats['Freq'] = row['freq']
            feats['num_edits'] = row['mismatches']
            feats['EditEff'] = row['Edit_Efficacy']
            all_features.append(feats)
            matched_gRNAS_out = pd.concat([matched_gRNAS_out,pd.DataFrame.from_dict(feats,orient='index').transpose()])

    matched_gRNAS_out = matched_gRNAS_out.fillna(0)
    matched_gRNAS_out.to_csv(output_csv, index=False)
    print(f"[✔] Features saved to {output_csv}")