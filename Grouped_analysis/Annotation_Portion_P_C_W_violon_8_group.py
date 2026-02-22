import os
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# === File paths ===
Bacteria_file_path = r'C:\Genomics\results\Genes_chardaff_parameters\Bacteria_2025-11-11_17-47-44.csv'
Archaea_file_path = r'C:\Genomics\results\Genes_chardaff_parameters\Archaea_2025-11-16_20-58-46.csv'
Fungi_file_path = r'C:\Genomics\results\Genes_chardaff_parameters\Fungisresults_2025-09-17_00-53-54.csv'
Arthropoda_file_path = r'C:\Genomics\results\Genes_chardaff_parameters\Arthropoda_500m.csv'
Chordata_file_path = r'C:\Genomics\results\Genes_chardaff_parameters\Chordata.csv'
Viride_file_path = r'C:\Genomics\results\Genes_chardaff_parameters\Viride_on_saveing_time_resuts.csv'
Mito_file_path = r'C:\Genomics\results\Genes_chardaff_parameters\Mito_2025-11-18_19-16-56.csv'
Chloro_file_path = r'C:\Genomics\results\Genes_chardaff_parameters\Chlroplast_2025-11-18_22-41-13.csv'

# === Load CSVs ===
Bacteria_df = pd.read_csv(Bacteria_file_path, index_col=0)
print(Bacteria_df.head().to_string())
Archaea_df = pd.read_csv(Archaea_file_path, index_col=0)
Fungi_df = pd.read_csv(Fungi_file_path, index_col=0)
Arthropoda_df = pd.read_csv(Arthropoda_file_path, index_col=0)
Chordata_df = pd.read_csv(Chordata_file_path, index_col=0)
Viride_df = pd.read_csv(Viride_file_path, index_col=0)
Mito_df = pd.read_csv(Mito_file_path, index_col=0)
Chloro_df = pd.read_csv(Chloro_file_path, index_col=0)

df = Mito_df


# === Column names ===
Wh_ano = df['Wh_ano']
P_anno = df['Int_P%'] + df['UTR_P%'] + df['CDS_P%'] + df['Share_P%']
C_anno = df['Int_C%'] + df['UTR_C%'] + df['CDS_C%'] + df['Share_C%']

Wh_ano = P_anno+C_anno

df = pd.DataFrame({'P_anno': P_anno, 'C_anno': C_anno,'Wh_anno':Wh_ano})
print(df.to_string())

# === STEP 4: Show statistics ===
print(f"\n--- {'Bacteria'} ---")
print(df.describe())

# === STEP 5: Plot ===
sns.set(style="whitegrid", font_scale=1.2)
fig, ax = plt.subplots(figsize=(10, 6))   # <-- single axis, not 4 subplots
#fig.suptitle(f"Nucleotide Composition Distributions - {'Bacteria'}", fontsize=14, weight="bold")

# Violin plot (all 4 categories in one figure)
sns.violinplot(
    data=df[['P_anno','C_anno','Wh_anno']],   # removed Wh_ano, keep only nucleotides
    ax=ax,
    palette=sns.light_palette("turquoise", n_colors=3),
    inner=None, linewidth=1.2, alpha=0.9
)

# Box overlay
sns.boxplot(
    data=df[['P_anno','C_anno','Wh_anno']],
    ax=ax,
    width=0.2,
    boxprops={'facecolor':'none', 'edgecolor':'black', 'linewidth':1.2},
    showcaps=False, showfliers=False,
    whiskerprops={'linewidth':0}
)

#ax.set_title("Nucleotide Composition (Violin + Box)")
#ax.set_ylabel("Percentage (%)")
ax.set_xticklabels(["P", "C", "Wh"],
                   fontsize=25, fontweight="bold")



# Calculate mean and std
means = df.mean()
stds  = df.std()

# Annotate mean values above each violin
for i, val in enumerate(means):
    ax.text(i, ax.get_ylim()[1]-5, f"μ:{val:.1f}", ha='center', va='bottom',
            color='black', fontsize=17, fontweight='bold',
            bbox=dict(facecolor='lightgray', alpha=0.6, boxstyle='round,pad=0.25'))

# Annotate std values just below the x-axis
for i, val in enumerate(stds):
    ax.text(i, ax.get_ylim()[0] + 2, f"σ:±{val:.1f}", ha='center', va='bottom',
            color='black', fontsize=17, fontweight='bold',
            bbox=dict(facecolor='lightgray', alpha=0.6, boxstyle='round,pad=0.25'))


plt.tight_layout()
#plt.savefig('Chloro_N_usage_violon.svg', dpi=900, format='svg',transparent=True)
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# === File paths ===
file_paths = {
    "Bacteria": r'C:\Genomics\results\Genes_chardaff_parameters\Bacteria_2025-11-11_17-47-44.csv',
    "Archaea": r'C:\Genomics\results\Genes_chardaff_parameters\Archaea_2025-11-16_20-58-46.csv',
    "Fungi": r'C:\Genomics\results\Genes_chardaff_parameters\Fungisresults_2025-09-17_00-53-54.csv',
    "Arthropoda": r'C:\Genomics\results\Genes_chardaff_parameters\Arthropoda_500m.csv',
    "Chordata": r'C:\Genomics\results\Genes_chardaff_parameters\Chordata.csv',
    "Viride": r'C:\Genomics\results\Genes_chardaff_parameters\Viride_on_saveing_time_resuts.csv',
    "Mito": r'C:\Genomics\results\Genes_chardaff_parameters\Mito_2025-11-18_19-16-56.csv',
    "Chloro": r'C:\Genomics\results\Genes_chardaff_parameters\Chlroplast_2025-11-18_22-41-13.csv'
}

# === Helper function to build dataframe ===
def build_df(df):
    P_anno = df['Int_P%'] + df['UTR_P%'] + df['CDS_P%'] + df['Share_P%']
    C_anno = df['Int_C%'] + df['UTR_C%'] + df['CDS_C%'] + df['Share_C%']
    Wh_anno = P_anno + C_anno
    return pd.DataFrame({'P_anno': P_anno, 'C_anno': C_anno, 'Wh_anno': Wh_anno})

# === Plotting ===
sns.set(style="whitegrid", font_scale=1.2)
fig, axes = plt.subplots(1, 8, figsize=(14, 3), sharey=True)  # 1 row, 8 columns

for ax, (name, path) in zip(axes, file_paths.items()):
    df_raw = pd.read_csv(path, index_col=0)
    df = build_df(df_raw)

    # Violin plot
    sns.violinplot(
        data=df[['P_anno','C_anno','Wh_anno']],
        ax=ax,
        palette=sns.light_palette("turquoise", n_colors=3),
        inner=None, linewidth=1.2, alpha=0.9
    )

    # Box overlay
    sns.boxplot(
        data=df[['P_anno','C_anno','Wh_anno']],
        ax=ax,
        width=0.2,
        boxprops={'facecolor':'none', 'edgecolor':'black', 'linewidth':1.2},
        showcaps=False, showfliers=False,
        whiskerprops={'linewidth':0}
    )

    # X labels
    ax.set_xticklabels(["P", "C", "Wh"], fontsize=12, fontweight="bold")
    ax.set_title(name, fontsize=14, weight="bold")

    # Annotate mean/std
    means = df.mean()
    stds  = df.std()
    for i, val in enumerate(means):
        ax.text(i, ax.get_ylim()[1]-5, f"μ:{val:.1f}", ha='center', va='bottom',
                color='black', fontsize=8, fontweight='bold',
                bbox=dict(facecolor='lightgray', alpha=0.6, boxstyle='round,pad=0.25'))
    for i, val in enumerate(stds):
        ax.text(i, 0, f"σ:±{val:.1f}", ha='center', va='bottom',
                color='black', fontsize=8, fontweight='bold',
                bbox=dict(facecolor='lightgray', alpha=0.6, boxstyle='round,pad=0.25'))

plt.tight_layout()
plt.savefig('grouped_whole_portion_and_strand_portion.svg',format='svg',dpi=900)
plt.show()

import matplotlib
print(matplotlib.get_backend())







