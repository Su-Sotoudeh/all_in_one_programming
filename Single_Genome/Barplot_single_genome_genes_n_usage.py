
import os
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import Counter
from Bio.Data import CodonTable

# === FONT SIZE CONFIGURATION (MASTER CONTROL) ===
TITLE_SIZE    = 20   # subplot title (folder name)
LABEL_SIZE    = 18   # axis labels
TICK_SIZE     = 18   # tick labels
ANNOT_SIZE    = 20   # % numbers above bars
SUPTITLE_SIZE = 16   # if you add a big title later
GENERAL_FONTSIZE = 18 # default body font

# Global settings (affects all plots)
plt.rcParams.update({
    "font.family": "serif",
    "font.size": GENERAL_FONTSIZE,
    "axes.titlesize": TITLE_SIZE,
    "axes.labelsize": LABEL_SIZE,
    "xtick.labelsize": TICK_SIZE,
    "ytick.labelsize": TICK_SIZE,
    "figure.titlesize": SUPTITLE_SIZE,
    "axes.linewidth": 0.8
})

sns.set(style="whitegrid")

# === PATH TO DATA ===
last_shered_path = r"C:\Genomics\Genomes\Bacteria_2025_2_1\bacteria_reference_complet_genomes\ncbi_dataset\data\\"
last_shered_path = r"C:\Genomics\Genomes\Penicilium\ncbi_dataset\data"
last_shered_path = r"C:\Genomics\Genomes\\fungi2"
last_shered_path = r"C:\Genomics\Genomes\Drosophila"
last_shered_path = r"C:\Genomics\Genomes\Human"
last_shered_path = r"C:\Genomics\Genomes\Arabidobsis"
last_shered_path = r"C:\Genomics\Genomes\Human_mito"
last_shered_path = r"C:\Genomics\Genomes\Arabidobsis_Chorlo"


list_of_bacteria_folders = os.listdir(last_shered_path)

for folder in list_of_bacteria_folders:
    if not folder.startswith("G"):
        continue
    #if folder.startswith("GCF_000007025.1"):
    files = os.listdir(os.path.join(last_shered_path, folder))
    fna = [f for f in files if (f.endswith(".fna") or f.endswith('.fasta')) and not f.startswith('cds')][0]
    gff = [f for f in files if f.endswith(".gff") or f.endswith('.gff3')][0]


    fna_file = os.path.join(last_shered_path, folder, fna)
    gff_file = os.path.join(last_shered_path, folder, gff)

    # === STEP 1: Read genome sequence ===
    genome = {}
    for record in SeqIO.parse(fna_file, "fasta"):
        genome[record.id] = record.seq.upper()

    # === STEP 2: Extract CDSs ===
    genes = []
    with open(gff_file) as gff_handle:
        for line in gff_handle:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            seqid, source, feature, start, end, score, strand, phase, attributes = parts

            if feature.lower() == 'cds':
                start, end = int(start), int(end)
                if seqid not in genome:
                    continue
                seq = genome[seqid][start-1:end]
                if strand == "-":
                    seq = seq.reverse_complement()
                genes.append(str(seq))

    # === STEP 3: Compute nucleotide composition ===
    genes_seq = ''.join(genes)
    len_genes = len(genes_seq)

    a = genes_seq.count("A") / len_genes * 100
    t = genes_seq.count("T") / len_genes * 100
    c = genes_seq.count("C") / len_genes * 100
    g = genes_seq.count("G") / len_genes * 100

    values = [a, t, c, g]
    labels = ['A', 'T', 'C', 'G']

    # Colorblind-friendly palette
    colors = ['#4daf4a', '#377eb8', '#ff7f00', '#b22222']

    # === PLOT ===
    fig, ax = plt.subplots(figsize=(6, 4.2))

    bars = ax.bar(labels, values, color=colors, edgecolor="black", linewidth=0.8)

    # % labels above bars
    for bar, v in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2,
                v + 0.5,
                f"{v:.2f}%",
                ha='center', va='bottom', fontsize=ANNOT_SIZE)

    #ax.set_ylabel("Proportion (%)", fontsize=LABEL_SIZE)
    ax.set_ylim(0, max(values) + 5)

    # Nature-style spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Horizontal grid only
    ax.yaxis.grid(True, linestyle="--", linewidth=0.4, alpha=0.6)
    ax.xaxis.grid(False)
    ax.set_axisbelow(True)
    ax.set_xticklabels(labels, fontsize=24)
    # Folder name as title
    #ax.set_title(folder, fontsize=TITLE_SIZE)

    plt.tight_layout()
    plt.savefig('Arabidobsis_Chorlo'+'_'+fna.split(".")[0] + ".svg", dpi=900,transparent=True)
    plt.show()








