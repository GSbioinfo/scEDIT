from Bio import SeqIO
import pandas as pd
import re
from collections import defaultdict

def reverse_complement(seq):
    return str(seq.translate(str.maketrans("ACGT", "TGCA"))[::-1])

def get_base_edit_window(seq, start=4, end=8):
    return seq[start-1:end]

def find_ng_pams(sequence):
    sequence = sequence.upper()
    valid_first_bases = {'A', 'T', 'C', 'G'}
    matches = []

    for i in range(len(sequence) - 1):
        dinuc = sequence[i:i+2]
        if dinuc[0] in valid_first_bases and dinuc[1] == 'G':
            matches.append((i, dinuc))

    return matches

def get_context_sequence(seq, pam_start, strand, total_length=30):
    """
    Extract 30nt context centered on gRNA + PAM. Adjusted for strand.
    """
    if strand == '+':
        center_start = max(pam_start - 20 - (total_length - 23)//2, 0)
        return seq[center_start:center_start + total_length]
    else:
        center_start = pam_start - 3 - (total_length - 23)//2
        center_start = max(center_start, 0)
        seq_segment = seq[center_start:center_start + total_length]
        return reverse_complement(seq_segment)

def find_gRNA_hits(seq1, gRNA, pam='NG', edit_window=(4, 8)):
    hits = []
    pam_regex = pam.replace('N', '[ACGT]')
    seq = seq1
    gRNARC= reverse_complement(gRNA)
    # Forward
    for match in re.finditer(pam_regex, seq):
        pam_start = match.start()
        if pam_start >= 20:
            candidate = seq[pam_start-20:pam_start]
            if candidate == gRNA:
                context = get_context_sequence(seq, pam_start, '+')
                hits.append(('+', pam_start, candidate, match.group(),
                             get_base_edit_window(candidate, *edit_window), context))

    # Reverse
    rc_seq = reverse_complement(seq1)
    for match in re.finditer(pam_regex, rc_seq):
        pam_start = match.start()
        if pam_start >= 20:
            candidate = rc_seq[pam_start-20:pam_start]
            if candidate == gRNA:
                orig_pos = len(seq) - pam_start - 2
                context = get_context_sequence(seq, orig_pos, '-')
                hits.append(('-', orig_pos, candidate, match.group(),
                             get_base_edit_window(candidate, *edit_window), context))
    
    return hits

def find_missing_gRNA(seq1, gRNA, pam='NG', edit_window=(4, 8)):
    hits = []
    pam_regex = pam.replace('N', '[ACGT]')
    seq = seq1
    gRNARC= reverse_complement(gRNA)
    rc_seq = reverse_complement(seq1)
    if(gRNA in seq):
        seq_site= seq.find(gRNA)
        seq_site_end = seq_site+len(gRNA)
        re.findall(pam_regex,seq[seq_site-5:seq_site_end+5])
    rc_seq_site = rc_seq.find(gRNA)
    seq_rcRNA_site= seq.find(gRNARC)
    rc_seq_rcRNA_site = rc_seq.find(gRNARC)
    for pos, pamx in find_ng_pams(seq):
        pam_start = pos
        if pam_start >= 20:
            candidate = seq[pam_start-20:pam_start]
            if candidate == gRNA:
                context = get_context_sequence(seq, pam_start, '+')
                hits.append(('+', pam_start, candidate, pamx,
                             get_base_edit_window(candidate, *edit_window), context))

    # Reverse
    rc_seq = reverse_complement(seq1)
    for pos, pamx in find_ng_pams(rc_seq):
        pam_start = pos
        if pam_start >= 20:
            candidate = rc_seq[pam_start-20:pam_start]
            if candidate == gRNARC:
                orig_pos = len(seq) - pam_start - 2
                context = get_context_sequence(seq, orig_pos, '-')
                hits.append(('-', orig_pos, candidate, pamx,
                             get_base_edit_window(candidate, *edit_window), context))
    
    return hits

def find_ng_pams_both_strands(sequence, gRNA, edit_window=(4, 8)):
    hits = []
    guide_length= len(gRNA)
    sequence = sequence.upper()
    rev_sequence = reverse_complement(sequence)
    gRNARC = reverse_complement(gRNA)
    valid_first_bases = {'A', 'T', 'C', 'G'}

    # Forward strand (+)
    for i in range(len(sequence) - 1):
        pam = sequence[i:i+2]
        if pam[0] in valid_first_bases and pam[1] == 'G':
            guide_start = max(i - guide_length, 0)
            guide_seq = sequence[guide_start:i] if i - guide_length >= 0 else "N/A"
            if(guide_seq == gRNA):
                context = get_context_sequence(sequence, i, '+')
                hits.append(('+', i, guide_seq, pam,
                             get_base_edit_window(guide_seq, *edit_window), context))


    # Reverse strand (−)
    for i in range(len(rev_sequence) - 1):
        pam = rev_sequence[i:i+2]
        if pam[0] in valid_first_bases and pam[1] == 'G':
            # Position in original sequence (reverse of reverse)
            orig_pos = len(sequence) - i - 2
            guide_end = orig_pos + 2 + guide_length
            guide_seq = sequence[orig_pos+2 : guide_end] if guide_end <= len(sequence) else "N/A"
            if(guide_seq == gRNA):
                context = get_context_sequence(sequence, orig_pos, '-')
                hits.append(('-', orig_pos, guide_seq, pam,
                             get_base_edit_window(guide_seq, *edit_window), context))

    return hits
    
def scan_fasta_all_guides(fasta_path, gRNA_file, pam='NG', edit_window=(4, 8), output_csv='output.csv'):
    gRNAs = set()
    #with open(gRNA_file) as f:
    #    for line in f:
    #        seq = line.strip().upper()
    #        if len(seq) == 20 and all(c in "ACGT" for c in seq):
    #            gRNAs.add(seq)
    gRNA_df= pd.read_csv(gRNA_file)
    for ro in gRNA_df.itertuples():
        seq = ro.gRNASeq.upper()
        if len(seq) == 20 and all(c in "ACGT" for c in seq):
            gRNAs.add(seq)
    results = []
    gRNA_counts = defaultdict(int)
    
    fastRecords = []
    fastaDict = {}
    jj=0
    for record in SeqIO.parse(fasta_path, "fasta"):
        fastRecords.append(record)
        fastaDict[record.id]=record
        jj=jj+1

    # Count hits
    #for record in SeqIO.parse(fasta_path, "fasta"):
    #    seq = str(record.seq).upper()
    #    for gRNA in gRNAs:
    #        gRNA_counts[gRNA] += len(find_gRNA_hits(seq, gRNA, pam, edit_window))
    for ro in gRNA_df.itertuples():
        seq = str(fastaDict[ro.AMPID].seq).upper()
        gRNA = ro.gRNASeq.upper()
        gRNA_counts[gRNA] += len(find_gRNA_hits(seq, gRNA, pam, edit_window))
    # Collect data
    #for record in SeqIO.parse(fasta_path, "fasta"):
    #    seq = str(record.seq).upper()
    #    for gRNA in gRNAs:
    #        hits = find_gRNA_hits(seq, gRNA, pam, edit_window)
    #        for strand, pam_pos, gseq, pam_seq, edit_seq, context_seq in hits:
    #            results.append({
    #                'chrom': record.id,
    #                'strand': strand,
    #                'gRNA': gseq,
    #                'pam': pam_seq,
    #                'pam_position': pam_pos,
    #                'edit_window_seq': edit_seq,
    #                'context_30nt': context_seq,
    #                'multiplicity': gRNA_counts[gRNA]
    #            })

    for ro in gRNA_df.itertuples():
        seq = str(fastaDict[ro.AMPID].seq).upper()
        gRNA = ro.gRNASeq.upper()
        hits = find_gRNA_hits(seq, gRNA, pam, edit_window)
        #hits = find_missing_gRNA(seq, gRNA, pam, edit_window)
        for strand, pam_pos, gseq, pam_seq, edit_seq, context_seq in hits:
            results.append({
                'chrom': ro.AMPID,
                'strand': strand,
                'gRNA': gseq,
                'pam': pam_seq,
                'pam_position': pam_pos,
                'edit_window_seq': edit_seq,
                'context_30nt': context_seq,
                'gRNA_rc':reverse_complement(gseq)
            })

    df = pd.DataFrame(results)
    #df.to_csv(output_csv, index=False)
    mapped = set(df['gRNA'].to_list())
    orig = set(gRNA_df['gRNASeq'].to_list())
    unmapped = list(sorted(orig - mapped))

    print(f"[x] Saved {len(df)} hits to {output_csv}")
    print(len(unmapped))
    unmapped_df = gRNA_df.loc[gRNA_df['gRNASeq'].isin(unmapped)]
    print(unmapped_df)
    for ro in unmapped_df.itertuples():
        seq = str(fastaDict[ro.AMPID].seq).upper()
        gRNA = ro.gRNASeq.upper()
        #hits = find_missing_gRNA(seq, gRNA, pam, edit_window)
        hits = find_ng_pams_both_strands(seq, gRNA, edit_window)
        for strand, pam_pos, gseq, pam_seq, edit_seq, context_seq in hits:
            results.append({
                'chrom': ro.AMPID,
                'strand': strand,
                'gRNA': gseq,
                'pam': pam_seq,
                'pam_position': pam_pos,
                'edit_window_seq': edit_seq,
                'context_30nt': context_seq,
                'gRNA_rc':reverse_complement(gseq)
            })
    newdf = pd.DataFrame(results)
    print(f"[x] Saved {len(newdf)} hits to {output_csv}")
    newmapped = set(newdf['gRNA'].to_list())
    newunmapped = list(sorted(orig - newmapped))
    print(len(newunmapped))
    newunmapped_df = gRNA_df.loc[gRNA_df['gRNASeq'].isin(newunmapped)]
    print(newunmapped_df)
    for ro in unmapped_df.itertuples():
        seq = str(fastaDict[ro.AMPID].seq).upper()
        gRNA = reverse_complement(ro.gRNASeq.upper())
        #hits = find_missing_gRNA(seq, gRNA, pam, edit_window)
        hits = find_ng_pams_both_strands(seq, gRNA, edit_window)
        for strand, pam_pos, gseq, pam_seq, edit_seq, context_seq in hits:
            results.append({
                'chrom': ro.AMPID,
                'strand': strand,
                'gRNA': reverse_complement(gseq),
                'pam': pam_seq,
                'pam_position': pam_pos,
                'edit_window_seq': edit_seq,
                'context_30nt': reverse_complement(context_seq),
                'gRNA_rc':gseq
            })
    newdf = pd.DataFrame(results)
    

    newdf["gRNAinContext"] = newdf.apply(lambda row: row["gRNA"].upper() not in row["context_30nt"].upper(), axis=1)
    mask = newdf.duplicated('gRNA', keep=False) & (newdf['gRNAinContext'])
    newdf_filtered = newdf[~mask]
    #newdf_filtered.loc[~newdf_filtered["gRNAinContext"],"context_30nt"]= reverse_complement(newdf_filtered["context_30nt"])
    newdf_filtered['context_30nt'] = newdf_filtered.apply(lambda row: reverse_complement(row["context_30nt"]) if row['gRNAinContext'] else row["context_30nt"], axis=1 )
    newdf_filtered["gRNAinContext"] = newdf_filtered.apply(lambda row: row["gRNA"].upper() not in row["context_30nt"].upper(), axis=1)
    
    final_filtered=newdf_filtered[newdf_filtered['context_30nt'].apply(lambda ro: len(ro)==30)].copy(deep=True)
    
    print(final_filtered['gRNAinContext'].value_counts())
    final_filtered.to_csv(output_csv, index=False)
    print(f"[✔] Saved {len( final_filtered)} hits to {output_csv}")
# ---------------------------
# Example usage
# ---------------------------
if __name__ == '__main__':
    scan_fasta_all_guides(
        fasta_path="/home/gajendra/Dropbox/Tapestri_workflow/PRJNA889818/temp_dir/Context_AMP_fullseq_GS.fasta" ,         # Replace with your FASTA
        gRNA_file='/home/gajendra/Dropbox/Tapestri_workflow/PRJNA889818/HBF1_gRNASeq_chrom_table.csv',              # TXT file, one 20-nt gRNA per line
        pam='NG',                           # evoFERNY PAM
        edit_window=(4, 8),                 # Adjust for base editor
        output_csv='HBF1_aiml_base_editor_targets.csv'
    )


'''print(f"[x] Saved {len(newdf)} hits to {output_csv}")
    newmapped = set(newdf['gRNA'].to_list())
    newunmapped = list(sorted(orig - newmapped))
    print(len(newunmapped))
    newunmapped_df = gRNA_df.loc[gRNA_df['gRNASeq'].isin(newunmapped)]
    print(newunmapped_df)
    
    for ro in unmapped_df.itertuples():
        seq = reverse_complement(str(fastaDict[ro.AMPID].seq).upper())
        gRNA = ro.gRNASeq.upper()
        #hits = find_missing_gRNA(seq, gRNA, pam, edit_window)
        hits = find_ng_pams_both_strands(seq, gRNA, edit_window)
        for strand, pam_pos, gseq, pam_seq, edit_seq, context_seq in hits:
            results.append({
                'chrom': ro.AMPID,
                'strand': strand,
                'gRNA': gseq,
                'pam': pam_seq,
                'pam_position': pam_pos,
                'edit_window_seq': edit_seq,
                'context_30nt': context_seq,
            })
    newdf = pd.DataFrame(results)
    print(f"[x] Saved {len(newdf)} hits to {output_csv}")
    newmapped = set(newdf['gRNA'].to_list())
    newunmapped = list(sorted(orig - newmapped))
    print(len(newunmapped))
    newunmapped_df = gRNA_df.loc[gRNA_df['gRNASeq'].isin(newunmapped)]
    print(newunmapped_df)'''