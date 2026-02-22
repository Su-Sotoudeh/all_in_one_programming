import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO

# File paths
sequence_file_path = r"/media/su/Bioinf/Genomes_2025_3/Genomes_2025_3_7/Chlroplast/ncbi_dataset/data/genome.fna"
df = pd.read_csv(r"/media/su/Bioinf/Genomes_2025_3/Genomes_2025_3_7/Chlroplast/ncbi_dataset/data/organelle_gene_report.csv")#, nrows=10)


# Define Okabe-Ito color palette (colorblind-friendly)
okabe_ito = ['#9467bd', '#2ca02c', '#1f77b4', '#ff7f0e', '#d62728',
             '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

# Lists to store results
P_portion_s = []
Aps = []
Tps = []
Cps = []
Gps = []
M_portion_s = []
Ams = []
Tms = []
Cms = []
Gms = []

# Process sequences
counter = 0
for record in SeqIO.parse(sequence_file_path, "fasta"):
    counter += 1
    if counter == 10:
        break
    print(f"Processing record {counter}")
    genes_sequence_plus = ''
    genes_sequence_minis = ''
    id = record.id
    sequence = record.seq

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
            continue
        total_length = len(genes_sequence_plus + genes_sequence_minis)
        P_portion = len(genes_sequence_plus) / total_length * 100 if total_length > 0 else 0
        Ap = genes_sequence_plus.count('A') / P_length * 100 if P_length > 0 else 0
        Tp = genes_sequence_plus.count('T') / P_length * 100 if P_length > 0 else 0
        Cp = genes_sequence_plus.count('C') / P_length * 100 if P_length > 0 else 0
        Gp = genes_sequence_plus.count('G') / P_length * 100 if P_length > 0 else 0
        M_portion = len(genes_sequence_minis) / total_length * 100 if total_length > 0 else 0
        Am = genes_sequence_minis.count('A') / M_length * 100 if M_length > 0 else 0
        Tm = genes_sequence_minis.count('T') / M_length * 100 if M_length > 0 else 0
        Cm = genes_sequence_minis.count('C') / M_length * 100 if M_length > 0 else 0
        Gm = genes_sequence_minis.count('G') / M_length * 100 if M_length > 0 else 0

        P_portion_s.append(P_portion)
        Aps.append(Ap)
        Tps.append(Tp)
        Cps.append(Cp)
        Gps.append(Gp)
        M_portion_s.append(M_portion)
        Ams.append(Am)
        Tms.append(Tm)
        Cms.append(Cm)
        Gms.append(Gm)

# Create a DataFrame for plotting
data = {
    'P_portion': P_portion_s,
    'A_plus': Aps,
    'T_plus': Tps,
    'C_plus': Cps,
    'G_plus': Gps,
    'M_portion': M_portion_s,
    'A_minus': Ams,
    'T_minus': Tms,
    'C_minus': Cms,
    'G_minus': Gms
}
plot_df = pd.DataFrame(data)
print(plot_df)
# Set style for better aesthetics
sns.set(style="whitegrid", font_scale=1.2)

# Create figure and axis
plt.figure(figsize=(15, 8))

# Create violin plot
ax = sns.violinplot(data=plot_df, palette=okabe_ito, inner="quartile", cut=0)

# Customize plot
plt.title("Distribution of Sequence Portions and Nucleotide Percentages", fontsize=16, pad=20)
plt.xlabel("Metrics", fontsize=12)
plt.ylabel("Percentage (%)", fontsize=12)

# Rotate x-axis labels for better readability
plt.xticks(rotation=45, ha='right')

# Calculate and display mean and std above each violin
for i, column in enumerate(plot_df.columns):
    mean_val = plot_df[column].mean()
    std_val = plot_df[column].std()
    max_val = plot_df[column].max()
    text = f'Mean: {mean_val:.2f}\nStd: {std_val:.2f}'
    ax.text(i, max_val + 5, text, ha='center', va='bottom', fontsize=10, color='black')

# Adjust layout to prevent label cutoff
plt.tight_layout()

# Save and show plot
plt.savefig("violin_plot.png", dpi=300, bbox_inches='tight')
plt.show()