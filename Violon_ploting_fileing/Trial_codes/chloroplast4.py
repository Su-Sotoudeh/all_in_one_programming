import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from Bio import SeqIO
import multiprocessing as mp
from functools import partial
import matplotlib.font_manager as fm

# Define expanded Okabe-Ito color palette (10 colorblind-friendly colors)
okabe_ito = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
             '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

# Function to process a single record (for parallelization)
def process_record(record, df):
    genes_sequence_plus = ''
    genes_sequence_minis = ''
    id = record.id
    sequence = record.seq
    result_p = {'P_portion': None, 'Ap': None, 'Tp': None, 'Cp': None, 'Gp': None}
    result_m = {'M_portion': None, 'Am': None, 'Tm': None, 'Cm': None, 'Gm': None}

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
            print(f"Skipping record {id}: No plus or minus strand genes")

            return None
        total_length = len(genes_sequence_plus) + len(genes_sequence_minis)
        P_portion = float(len(genes_sequence_plus) / total_length * 100 if total_length > 0 else 0)
        Ap = float(genes_sequence_plus.count('A') / P_length * 100 if P_length > 0 else 0)
        Tp = float(genes_sequence_plus.count('T') / P_length * 100 if P_length > 0 else 0)
        Cp = float(genes_sequence_plus.count('C') / P_length * 100 if P_length > 0 else 0)
        Gp = float(genes_sequence_plus.count('G') / P_length * 100 if P_length > 0 else 0)
        M_portion = float(len(genes_sequence_minis) / total_length * 100 if total_length > 0 else 0)
        Am = float(genes_sequence_minis.count('A') / M_length * 100 if M_length > 0 else 0)
        Tm = float(genes_sequence_minis.count('T') / M_length * 100 if M_length > 0 else 0)
        Cm = float(genes_sequence_minis.count('C') / M_length * 100 if M_length > 0 else 0)
        Gm = float(genes_sequence_minis.count('G') / M_length * 100 if M_length > 0 else 0)
        result_p = {'P_portion': P_portion, 'Ap': Ap, 'Tp': Tp, 'Cp': Cp, 'Gp': Gp}
        result_m = {'M_portion': M_portion, 'Am': Am, 'Tm': Tm, 'Cm': Cm, 'Gm': Gm}
        # Verify all values are numeric
        if any(v is None or not isinstance(v, (int, float)) for v in result_p.values()) or \
                any(v is None or not isinstance(v, (int, float)) for v in result_m.values()):
            print(f"Skipping record {id}: Non-numeric values in results")
            return None
    return result_p, result_m

# Data processing
sequence_file_path = r"/media/su/Bioinf/Genomes_2025_3/Genomes_2025_3_7/Chlroplast/ncbi_dataset/data/genome.fna"
df = pd.read_csv(r"/media/su/Bioinf/Genomes_2025_3/Genomes_2025_3_7/Chlroplast/ncbi_dataset/data/organelle_gene_report.csv")
print("CSV head:")
print(df.head())

# Parallel processing (limit to 10 valid records)
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

print("Starting parallel processing...")
with mp.Pool(processes=mp.cpu_count()) as pool:
    records = list(SeqIO.parse(sequence_file_path, "fasta"))
    results = pool.map(partial(process_record, df=df), records)

# Collect results (up to 10 valid records)
counter = 0
for result in results:
    if result is not None and counter < 10:
        print(f"Processing record {result[0].get('P_portion', 'N/A')}, {result[1].get('M_portion', 'N/A')}")
        # Check for valid numeric values
        if (isinstance(result[0]['P_portion'], (int, float)) and np.isfinite(result[0]['P_portion']) and
                isinstance(result[0]['Ap'], (int, float)) and np.isfinite(result[0]['Ap']) and
                isinstance(result[0]['Tp'], (int, float)) and np.isfinite(result[0]['Tp']) and
                isinstance(result[0]['Cp'], (int, float)) and np.isfinite(result[0]['Cp']) and
                isinstance(result[0]['Gp'], (int, float)) and np.isfinite(result[0]['Gp']) and
                isinstance(result[1]['M_portion'], (int, float)) and np.isfinite(result[1]['M_portion']) and
                isinstance(result[1]['Am'], (int, float)) and np.isfinite(result[1]['Am']) and
                isinstance(result[1]['Tm'], (int, float)) and np.isfinite(result[1]['Tm']) and
                isinstance(result[1]['Cm'], (int, float)) and np.isfinite(result[1]['Cm']) and
                isinstance(result[1]['Gm'], (int, float)) and np.isfinite(result[1]['Gm'])):
            counter += 1
            print(f"Processed record {counter}")
            P_portion_s.append(result[0]['P_portion'])
            Aps.append(result[0]['Ap'])
            Tps.append(result[0]['Tp'])
            Cps.append(result[0]['Cp'])
            Gps.append(result[0]['Gp'])
            M_portion_s.append(result[1]['M_portion'])
            Ams.append(result[1]['Am'])
            Tms.append(result[1]['Tm'])
            Cms.append(result[1]['Cm'])
            Gms.append(result[1]['Gm'])
        else:
            print(f"Skipping record: Invalid values - {result}")
    if counter >= 10:
        break

# Verify data
if not P_portion_s:
    raise ValueError("No valid records processed. Check accession IDs in CSV and FASTA file or remove nrows limit.")

# Create DataFrame for 10 groups
print('Chloroplast',len(Aps), len(Tps), len(Cps), len(Gms),len(Ams), len(Tms), len(Cms))
plot_data = pd.DataFrame({
    "value": (P_portion_s + Aps + Tps + Cps + Gps +
              M_portion_s + Ams + Tms + Cms + Gms),
    "group": (["P_portion"] * len(P_portion_s) +
              ["Ap"] * len(Aps) +
              ["Tp"] * len(Tps) +
              ["Cp"] * len(Cps) +
              ["Gp"] * len(Gps) +
              ["M_portion"] * len(M_portion_s) +
              ["Am"] * len(Ams) +
              ["Tm"] * len(Tms) +
              ["Cm"] * len(Cms) +
              ["Gm"] * len(Gms))
})

# Check for NaN/Inf in plot_data
if plot_data["value"].isna().any() or np.isinf(plot_data["value"]).any():
    raise ValueError("plot_data contains NaN or Inf values. Check data processing.")

# Check available fonts and select Arial or fallback
available_fonts = [f.name for f in fm.fontManager.ttflist]
font_family = 'Arial' if 'Arial' in available_fonts else 'DejaVu Sans'
print(f"Using font: {font_family}")

# Set custom style for Nature Communications
plt.style.use('default')
plt.rcParams.update({
    'font.family': font_family,
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
    'axes.grid': False,
})

# Create violin plot for 10 groups
fig, ax = plt.subplots(figsize=(10, 5))  # Wider figure for 10 groups
sns.violinplot(x="group", y="value", hue="group", data=plot_data, ax=ax, palette=okabe_ito, linewidth=0.8, inner="box", legend=False)
ax.set_title('Mitochondrial Genome Strand and Nucleotide Composition', fontsize=12)
ax.set_xlabel('Metric', labelpad=5)
ax.set_ylabel('Percentage (%)', labelpad=5)

# Add mean and SD annotations
stats = plot_data.groupby("group")["value"].agg(["mean", "std"]).reset_index()
for i, group in enumerate(plot_data["group"].unique()):
    m = stats.loc[stats["group"] == group, "mean"].values[0]
    s = stats.loc[stats["group"] == group, "std"].values[0]
    ymax = plot_data[plot_data["group"] == group]["value"].max()
    ax.text(i, ymax * 1.2, f"Mean={m:.2f}\nSD={s:.2f}",
            ha='center', va='bottom', fontsize=8, color="black")

# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Set y-axis limit with fallback for edge cases
max_value = plot_data["value"].max()
if np.isfinite(max_value):
    ax.set_ylim(0,100 )#max_value * 1.2)
else:
    ax.set_ylim(0, 100)  # Fallback to reasonable range

# Rotate x-axis labels for readability with 10 groups
plt.xticks(rotation=45, ha='right')

# Save figure
plt.tight_layout()
plt.savefig("mtgenome_violin_ten_groups.png", dpi=600, bbox_inches="tight", format='png')
plt.savefig("mtgenome_violin_ten_groups.tiff", dpi=600, bbox_inches="tight", format='tiff')
plt.savefig("mtgenome_violin_ten_groups.eps", bbox_inches="tight", format='eps')

plt.show()

# Print statistics
print("\nStatistics:")
print(stats)