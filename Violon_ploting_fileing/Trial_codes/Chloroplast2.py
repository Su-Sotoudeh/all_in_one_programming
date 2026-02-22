import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from Bio import SeqIO
import multiprocessing as mp
from functools import partial

# Define Okabe-Ito color palette (colorblind-friendly)
okabe_ito = ['#9467bd', '#2ca02c','#1f77b4', '#ff7f0e','#d62728', ]  # Blue, Orange, Green, Red, Purple

# Function to process a single record (for parallelization)
def process_record(record, df):
    genes_sequence_plus = ''
    genes_sequence_minis = ''
    id = record.id
    sequence = record.seq
    result_p = {'P_poerion': None, 'Ap': None, 'Tp': None, 'Cp': None, 'Gp': None}
    result_m = {'M_portion': None, 'Ac': None, 'Tc': None, 'Cc': None, 'Gc': None}

    if id.startswith("NC"):
        one_mt_df = df[df["accession"] == id]
        plus = 0
        minis = 0
        P_length = 0
        M_length = 0
        for _, row in one_mt_df.iterrows():
            start = row["start"] - 1
            end = row["end"]
            gene_sequence = sequence[start:end]
            length = end - start
            if row["orientation"] == "plus":
                genes_sequence_plus += gene_sequence
                plus += 1
                P_length += length
            else:
                genes_sequence_minis += gene_sequence
                minis += 1
                M_length += length
        if plus == 0 and minis == 0:
            return None
        P_portion = genes_sequence_plus / (genes_sequence_plus + genes_sequence_minis) * 100 if (genes_sequence_plus + genes_sequence_minis) > 0 else 0
        Ap = genes_sequence_plus.count('A') / P_length * 100 if P_length > 0 else 0
        Tp = genes_sequence_plus.count('T') / P_length * 100 if P_length > 0 else 0
        Cp = genes_sequence_plus.count('C') / P_length * 100 if P_length > 0 else 0
        Gp = genes_sequence_plus.count('G') / P_length * 100 if P_length > 0 else 0
        result_p = {'P_portion': P_portion, 'A': Ap, 'T': Tp, 'C': Cp, 'G': Gp}
        M_portion = genes_sequence_minis / (genes_sequence_plus + genes_sequence_minis) * 100 if (genes_sequence_plus + genes_sequence_minis) > 0 else 0
        Am = genes_sequence_plus.count('A') / P_length * 100 if P_length > 0 else 0
        Tm = genes_sequence_plus.count('T') / P_length * 100 if P_length > 0 else 0
        Cm = genes_sequence_plus.count('C') / P_length * 100 if P_length > 0 else 0
        Gm = genes_sequence_plus.count('G') / P_length * 100 if P_length > 0 else 0
        result_m = {'M_portion': M_portion, 'A': Am, 'T': Tm, 'C': Cm, 'G': Gm}

    return result_p,result_m

# Data processing
sequence_file_path = r"/media/su/Bioinf/Genomes_2025_3/Genomes_2025_3_7/Chlroplast/ncbi_dataset/data/genome.fna"
df = pd.read_csv(r"/media/su/Bioinf/Genomes_2025_3/Genomes_2025_3_7/Chlroplast/ncbi_dataset/data/organelle_gene_report.csv")

# Parallel processing
P2C_length_list = []
As = []
Ts = []
Cs = []
Gs = []

print("Starting parallel processing...")
with mp.Pool(processes=mp.cpu_count()) as pool:
    records = list(SeqIO.parse(sequence_file_path, "fasta"))
    results = pool.map(partial(process_record, df=df), records)

# Collect results
counter = 0
for result in results:
    if result is not None:
        counter += 1
        print(f"Processed record {counter}")
        P2C_length_list.append(result['P2C_length'])
        As.append(result['A'])
        Ts.append(result['T'])
        Cs.append(result['C'])
        Gs.append(result['G'])

# Create DataFrame
plot_data = pd.DataFrame({
    "value": P2C_length_list + As + Ts + Cs + Gs,
    "group": (["P/M"] * len(P2C_length_list) +
              ["A"] * len(As) +
              ["T"] * len(Ts) +
              ["C"] * len(Cs) +
              ["G"] * len(Gs))
})

# Set custom style for Nature Communications
plt.style.use('default')  # Reset to default for full control
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 10,
    'axes.labelsize': 10,
    'axes.titlesize': 12,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'figure.dpi': 300,
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'axes.facecolor': 'white',
    'figure.facecolor': 'white',
    'axes.grid': False,  # No grid for cleaner look
})

# Create violin plot
fig, ax = plt.subplots(figsize=(5, 4))
sns.violinplot(x="group", y="value", data=plot_data, ax=ax, palette=okabe_ito, linewidth=0.8, inner="box")
ax.set_title('Mitochondrial Genome Nucleotide Composition', fontsize=12)
ax.set_xlabel('Nucleotide/Strand', labelpad=5)
ax.set_ylabel('Percentage (%)', labelpad=5)

# Add mean and SD annotations
stats = plot_data.groupby("group")["value"].agg(["mean", "std"]).reset_index()
for i, group in enumerate(plot_data["group"].unique()):
    m = stats.loc[stats["group"] == group, "mean"].values[0]
    s = stats.loc[stats["group"] == group, "std"].values[0]
    ymax = plot_data[plot_data["group"] == group]["value"].max()
    ax.text(i, ymax * 1.05, f"Mean={m:.2f}\nSD={s:.2f}",
            ha='center', va='bottom', fontsize=8, color="black")

# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Set y-axis limit
ax.set_ylim(0,100)# max(plot_data["value"]) )

# Save figure
plt.tight_layout()
plt.savefig("mtgenome_violin.png", dpi=600, bbox_inches="tight", format='png')
plt.savefig("mtgenome_violin.tiff", dpi=600, bbox_inches="tight", format='tiff')
plt.savefig("mtgenome_violin.eps", bbox_inches="tight", format='eps')

plt.show()

# Print statistics
print("\nStatistics:")
print(stats)