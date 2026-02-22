import os
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import Counter
from Bio.Data import CodonTable

# === FONT SIZE CONFIGURATION ===
TITLE_SIZE    = 11   # subplot titles
LABEL_SIZE    = 10   # axis labels
TICK_SIZE     = 10   # tick labels
LEGEND_SIZE   = 10  # legend text
SUPTITLE_SIZE = 11   # figure title

# Apply globally
plt.rcParams['axes.titlesize']   = TITLE_SIZE
plt.rcParams['axes.labelsize']   = LABEL_SIZE
plt.rcParams['xtick.labelsize']  = TICK_SIZE
plt.rcParams['ytick.labelsize']  = TICK_SIZE
plt.rcParams['legend.fontsize']  = LEGEND_SIZE
plt.rcParams['figure.titlesize'] = SUPTITLE_SIZE

sns.set(style="whitegrid", font_scale=1.2)  # optional scaling

# === PATH TO DATA ===
last_shered_path = r"C:\Genomics\Genomes\Bacteria_2025_2_1\bacteria_reference_complet_genomes\ncbi_dataset\data\\"
last_shered_path = r"C:\Genomics\Genomes\Penicilium\ncbi_dataset\data"
last_shered_path = r"C:\Genomics\Genomes\fungi2"
#last_shered_path = r"C:\Genomics\Genomes\Drosophila"
#last_shered_path = r"C:\Genomics\Genomes\Human"
#last_shered_path = r"C:\Genomics\Genomes\Arabidobsis"
#last_shered_path = r"C:\Genomics\Genomes\Human_mito"
#last_shered_path = r"C:\Genomics\Genomes\Arabidobsis_Chorlo"
list_of_bacteria_folders = os.listdir(last_shered_path)
print(list_of_bacteria_folders)
for folder in list_of_bacteria_folders:
    if not folder.startswith("G"):
        continue

    files = os.listdir(os.path.join(last_shered_path, folder))
    fna = [f for f in files if (f.endswith(".fna") or f.endswith('.fasta')) and not f.startswith('cds')][0]
    gff = [f for f in files if f.endswith(".gff") or f.endswith('.gff3')][0]
    print(fna, gff)
    fna_file = os.path.join(last_shered_path, folder, fna)
    gff_file = os.path.join(last_shered_path, folder, gff)

    # === STEP 1: Read genome sequence ===
    genome = {}
    for record in SeqIO.parse(fna_file, "fasta"):
        genome[record.id] = record.seq.upper()

    # === STEP 2: Extract gene regions from GFF ===
    genes = []
    with open(gff_file) as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            seqid, source, feature, start, end, score, strand, phase, attributes = parts

            if feature.lower() == 'gene':
                start, end = int(start), int(end)
                if seqid not in genome:
                    continue
                seq = genome[seqid][start-1:end]
                if strand == "-":
                    seq = seq.reverse_complement()
                genes.append(str(seq))

    # === STEP 3: Calculate nucleotide composition per gene ===
    data = []
    for i, gene_seq in enumerate(genes, 1):
        gene_len = len(gene_seq)
        if gene_len == 0:
            continue
        a = gene_seq.count("A") / gene_len * 100
        t = gene_seq.count("T") / gene_len * 100
        c = gene_seq.count("C") / gene_len * 100
        g = gene_seq.count("G") / gene_len * 100
        at_ratio = a / t if t != 0 else np.nan
        gc_ratio = g / c if c != 0 else np.nan
        GC_count = (g+c)/(a+t+c+g)*100
        #Annotation_portion =
        data.append([i, a, t, c, g, at_ratio, gc_ratio, GC_count])

    df = pd.DataFrame(data, columns=["Gene", "A%", "T%", "C%", "G%", "A/T", "G/C", 'GC%'])

    # === STEP 4: Show statistics ===
    print(f"\n--- {folder} ---")
    print(df.describe())

    # === STEP 5: Plot ===
    fig, axes = plt.subplots(3, 1, figsize=(3, 7))
    #fig.suptitle(f"Nucleotide Composition Distributions - {folder}", fontsize=SUPTITLE_SIZE, weight="bold")

    # --- VIOLIN 1 ---
    ax = axes[0]
    sns.violinplot(data=df[["A%", "T%", "C%", "G%"]],
                   ax=ax,
                   palette=sns.light_palette("turquoise", n_colors=4),
                   inner=None, linewidth=1.2, alpha=0.9)
    sns.boxplot(data=df[["A%", "T%", "C%", "G%"]],
                ax=ax,
                width=0.2, boxprops={'facecolor':'none', 'edgecolor':'black', 'linewidth':1.2},
                showcaps=False, showfliers=False, whiskerprops={'linewidth':0})
    ax.set_title("Nucleotide Composition (Violin + Box)", fontsize=TITLE_SIZE)
    #ax.set_ylabel("Percentage (%)", fontsize=LABEL_SIZE)
    ax.set_ylim(0, df[["A%", "T%", "C%", "G%"]].max().max()*1.1)
    ax.grid(True)
    ax.tick_params(axis='both', labelsize=TICK_SIZE)

    means = df[["A%", "T%", "C%", "G%"]].mean()
    for i, val in enumerate(means):
        ax.text(i, 0.5, f"₰:{val:.1f}", ha='center', va='center',
                color='black', fontsize=8, fontweight='bold',
                bbox=dict(facecolor='white', alpha=0.6, boxstyle='round,pad=0.25'))

    # --- SCATTER 2 ---
    ax = axes[1]
    ax.set_title("Nucleotide Usage in Diffrence GC%", fontsize=TITLE_SIZE)
    ax.scatter(df['GC%'], df['A%'], s=0.5, color='green')
    ax.scatter(df['GC%'], df['T%'], s=0.5, color='blue')
    ax.scatter(df['GC%'], df['C%'], s=0.5, color='yellow')
    ax.scatter(df['GC%'], df['G%'], s=0.5, color='red')
    ax.set_xlabel("GC%", fontsize=LABEL_SIZE)
    #ax.set_ylabel("Percentage (%)", fontsize=LABEL_SIZE)
    ax.set_ylim(0, df[["A%", "T%", "C%", "G%"]].max().max()*1.1)
    ax.grid(True)
    ax.tick_params(axis='both', labelsize=TICK_SIZE)

    # --- VIOLIN 3 ---
    ax = axes[2]
    sns.violinplot(data=df[["A/T", "G/C"]],
                   ax=ax,
                   palette=sns.light_palette("teal", n_colors=2),
                   inner=None, linewidth=1.2, alpha=0.9)
    sns.boxplot(data=df[["A/T", "G/C"]],
                ax=ax,
                width=0.2, boxprops={'facecolor':'none', 'edgecolor':'black', 'linewidth':1.2},
                showcaps=False, showfliers=False, whiskerprops={'linewidth':0})
    ax.set_title("A/T and G/C Ratios (Violin + Box)", fontsize=TITLE_SIZE)
    #ax.set_ylabel("Ratio", fontsize=LABEL_SIZE)
    ax.set_ylim(0, 4)# df[["G/C",'A/T']].max().max()*1.1)
    ax.grid(True)
    ax.tick_params(axis='both', labelsize=TICK_SIZE)

    means = df[["A/T", "G/C"]].mean()
    for i, val in enumerate(means):
        ax.text(i, 0, f"₰{val:.2f}", ha='center', va='center',
                color='black', fontsize=8, fontweight='bold',
                bbox=dict(facecolor='white', alpha=0.6, boxstyle='round,pad=0.25'))




    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig("chloro_n_composition.png",transparent=True,dpi=900)
    plt.show()