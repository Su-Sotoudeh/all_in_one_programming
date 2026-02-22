import matplotlib
matplotlib.use("Agg")

import os
from Bio import SeqIO
import pandas as pd
import numpy as np
from itertools import product
from concurrent.futures import ProcessPoolExecutor

# === Paths ===
data_path = r"C:\Genomics\Genomes\Bacteria_2025_2_1\bacteria_reference_complet_genomes\ncbi_dataset\data\\"
list_of_bacteria_folders = os.listdir(data_path)
list_of_bacteria_folders = list_of_bacteria_folders

# === All 64 codons ===
bases = ['A','T','C','G']
codons = [''.join(p) for p in product(bases, repeat=3)]

def process_folder(folder):
    """Process one bacterial genome folder and return a row dict."""
    if not folder.startswith("G"):
        return None

    files = os.listdir(os.path.join(data_path, folder))
    if len(files) < 2:
        return None

    fna_files = [f for f in files if f.endswith(".fna")]
    gff_files = [f for f in files if f.endswith(".gff")]
    if not fna_files or not gff_files:
        return None

    fna_path = os.path.join(data_path, folder, fna_files[0])
    gff_path = os.path.join(data_path, folder, gff_files[0])

    # === Read genome ===
    genome = {}
    for record in SeqIO.parse(fna_path, "fasta"):
        genome[record.id] = record.seq.upper()

    # === Extract CDS sequences ===
    cds_list = []
    with open(gff_path) as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            seqid, source, feature, start, end, score, strand, phase, attributes = parts
            if feature.lower() != "cds":
                continue
            start, end = int(start), int(end)
            if seqid not in genome:
                continue
            seq = genome[seqid][start-1:end]
            if strand == "-":
                seq = seq.reverse_complement()
            cds_list.append(str(seq))

    # === Count codons ===
    codon_counts = dict.fromkeys(codons, 0)
    total_codons = 0
    for cds in cds_list:
        for i in range(0, len(cds)-2, 3):
            codon = cds[i:i+3]
            if codon in codon_counts:
                codon_counts[codon] += 1
                total_codons += 1

    # === Normalize codon frequencies ===
    codon_freq = {k: (v/total_codons)*100 if total_codons > 0 else 0 for k, v in codon_counts.items()}

    # === Compute genome GC% ===
    all_seq = ''.join(cds_list)
    g = all_seq.count('G')
    c = all_seq.count('C')
    a = all_seq.count('A')
    t = all_seq.count('T')
    total_nt = g + c + a + t
    gc_percent = (g + c) / total_nt * 100 if total_nt > 0 else np.nan

    # === Return row dict ===
    row = {'Organism': folder, 'GC%': gc_percent}
    row.update(codon_freq)
    return row

# === Run in parallel ===
results = []
with ProcessPoolExecutor() as executor:
    for i, row in enumerate(executor.map(process_folder, list_of_bacteria_folders), start=1):
        if row is not None:
            results.append(row)
            print(i, row['Organism'])

# === Build dataframe and save ===
codon_usage_df = pd.DataFrame(results, columns=['Organism', 'GC%'] + codons)
codon_usage_df.to_csv("codon_usage_df_6000.csv", index=False)
print(codon_usage_df.head().to_string())

os.system('paplay /usr/share/sounds/freedesktop/stereo/complete.oga')

