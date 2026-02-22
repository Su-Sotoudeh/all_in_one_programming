import os
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter, defaultdict

from collections import Counter
from Bio.Data import CodonTable
import matplotlib.pyplot as plt







# === PATH TO DATA ===
last_shered_path = r"C:\Genomics\Genomes\Bacteria_2025_2_1\bacteria_reference_complet_genomes\ncbi_dataset\data\\"
list_of_bacteria_folders = os.listdir(last_shered_path)

for folder in list_of_bacteria_folders:
    if not folder.startswith("G"):
        continue

    files = os.listdir(os.path.join(last_shered_path, folder))
    fna = [f for f in files if f.endswith(".fna")][0]
    gff = [f for f in files if f.endswith(".gff")][0]
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

            if feature.lower() == 'cds':
                start, end = int(start), int(end)
                if seqid not in genome:
                    continue
                seq = genome[seqid][start-1:end]
                if strand == "-":
                    seq = seq.reverse_complement()
                genes.append(str(seq))






    # === STEP 3: Calculate nucleotide composition per gene ===
    from itertools import product
    from collections import Counter
    import pandas as pd
    import numpy as np

    # -------------------------------
    # Your original data structure
    # -------------------------------
    data = []

    # Independent 64 codons (canonical order)
    bases = ["A", "T", "G", "C"]
    all_codons = [''.join(p) for p in product(bases, repeat=3)]
    codon_columns = [f"{codon}" for codon in all_codons]

    for i, gene_seq in enumerate(genes, 1):
        gene_seq = gene_seq.upper().replace("\n", "").replace(" ", "")
        gene_len = len(gene_seq)
        if gene_len == 0:
            continue

        # --- nucleotide composition ---
        a = gene_seq.count("A") / gene_len * 100
        t = gene_seq.count("T") / gene_len * 100
        c = gene_seq.count("C") / gene_len * 100
        g = gene_seq.count("G") / gene_len * 100
        L = gene_seq.count("A")+  gene_seq.count("T")+ gene_seq.count("C")+ gene_seq.count("G")
        at_ratio = a / t if t != 0 else np.nan
        gc_ratio = g / c if c != 0 else np.nan
        GC_count = (g + c) / (a + t + c + g) * 100

        # --- codon counts ---
        n_codons = gene_len // 3
        codons = [gene_seq[j*3:(j*3)+3] for j in range(n_codons)]
        codon_count = Counter(codons)

        # Build row
        row = [i, L,a, t, c, g, at_ratio, gc_ratio, GC_count]
        # Add all 64 codons in order, zero if missing
        row += [codon_count.get(codon, 0) for codon in all_codons]

        data.append(row)

    # --- Create DataFrame ---
    df = pd.DataFrame(
        data,
        columns=["Gene",'L', "A%", "T%", "C%", "G%", "A/T", "G/C", "GC%"] + codon_columns
    )

    print(df.head().to_string())








    # === STEP 4: Show statistics ===
    print(f"\n--- {folder} ---")
    print(df.describe())

    # === STEP 5: Plot ===
    sns.set(style="whitegrid", font_scale=1.2)
    fig, axes = plt.subplots(3, 1, figsize=(2.5,8))
    fig.suptitle(f"Nucleotide Composition Distributions - {folder}", fontsize=14, weight="bold")



    # ---------------- VIOLIN 1: A, T, C, G (with box inside) ----------------
    ax = axes[0]
    sns.violinplot(data=df[["A%", "T%", "C%", "G%"]],
                   ax=ax,
                   palette=sns.light_palette("turquoise", n_colors=4),
                   inner=None, linewidth=1.2, alpha=0.9)
    sns.boxplot(data=df[["A%", "T%", "C%", "G%"]],
                ax=ax,
                width=0.2, boxprops={'facecolor':'none', 'edgecolor':'black', 'linewidth':1.2},
                showcaps=False, showfliers=False, whiskerprops={'linewidth':0})
    ax.set_title("Nucleotide Composition (Violin + Box)")
    ax.set_ylabel("Percentage (%)")

    means = df[["A%", "T%", "C%", "G%"]].mean()
    min =   df[["A%", "T%", "C%", "G%"]].min()
    for i, val in enumerate(means):
        ax.text(i,0.5, f"₰:{val:.1f}", ha='center', va='center',
                color='black', fontsize=7, fontweight='bold',
                bbox=dict(facecolor='white', alpha=0.6, boxstyle='round,pad=0.25'))

    # ---------------- scatter 2: A,T, G,C  ----------------
    ax = axes[1]

    ax.scatter(df['GC%'], df['A%'],s = 0.5,color = 'green')
    ax.scatter(df['GC%'], df['T%'],s = 0.5,color = 'blue')
    ax.scatter(df['GC%'], df['C%'],s = 0.5,color = 'yellow')
    ax.scatter(df['GC%'], df['G%'],s = 0.5,color = 'red')




    # ---------------- VIOLIN 3: A/T, G/C (with box inside) ----------------
    ax = axes[2]
    sns.violinplot(data=df[["A/T", "G/C"]],
                   ax=ax,
                   palette=sns.light_palette("teal", n_colors=2),
                   inner=None, linewidth=1.2, alpha=0.9)
    sns.boxplot(data=df[["A/T", "G/C"]],
                ax=ax,
                width=0.2, boxprops={'facecolor':'none', 'edgecolor':'black', 'linewidth':1.2},
                showcaps=False, showfliers=False, whiskerprops={'linewidth':0})
    ax.set_title("A/T and G/C Ratios (Violin + Box)")
    ax.set_ylabel("Ratio")

    means = df[["A/T", "G/C"]].mean()
    for i, val in enumerate(means):
        ax.text(i, 0.5, f"₰{val:.2f}", ha='center', va='center',
                color='black', fontsize=7, fontweight='bold',
                bbox=dict(facecolor='white', alpha=0.6, boxstyle='round,pad=0.25'))






    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.show()
