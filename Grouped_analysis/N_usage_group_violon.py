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
Archaea_df = pd.read_csv(Archaea_file_path, index_col=0)
Fungi_df = pd.read_csv(Fungi_file_path, index_col=0)
Arthropoda_df = pd.read_csv(Arthropoda_file_path, index_col=0)
Chordata_df = pd.read_csv(Chordata_file_path, index_col=0)
Viride_df = pd.read_csv(Viride_file_path, index_col=0)
Mito_df = pd.read_csv(Mito_file_path, index_col=0)
Chloro_df = pd.read_csv(Chloro_file_path, index_col=0)

# === Column names ===
Wh_ano = 'Wh_ano'
A = 'CDS_P_A_ano'
T = 'CDS_P_T_ano'
C = 'CDS_P_C_ano'
G = 'CDS_P_G_ano'

df = Chloro_df

# === STEP 4: Show statistics ===
print(f"\n--- {'Bacteria'} ---")
print(df.describe())

# === STEP 5: Plot ===
sns.set(style="whitegrid", font_scale=1.2)
fig, ax = plt.subplots(figsize=(10, 6))   # <-- single axis, not 4 subplots
#fig.suptitle(f"Nucleotide Composition Distributions - {'Bacteria'}", fontsize=14, weight="bold")

# Violin plot (all 4 categories in one figure)
sns.violinplot(
    data=df[[A, T, C, G]],   # removed Wh_ano, keep only nucleotides
    ax=ax,
    palette=sns.light_palette("turquoise", n_colors=4),
    inner=None, linewidth=1.2, alpha=0.9
)

# Box overlay
sns.boxplot(
    data=df[[A, T, C, G]],
    ax=ax,
    width=0.2,
    boxprops={'facecolor':'none', 'edgecolor':'black', 'linewidth':1.2},
    showcaps=False, showfliers=False,
    whiskerprops={'linewidth':0}
)

#ax.set_title("Nucleotide Composition (Violin + Box)")
#ax.set_ylabel("Percentage (%)")
ax.set_xticklabels(["A", "T", "C", "G"],
                   fontsize=25, fontweight="bold")
'''# Annotate means
means = df[[A, T, C, G]].mean()
for i, val in enumerate(means):
    ax.text(i, val, f"mean:{val:.1f}", ha='center', va='center',
            color='black', fontsize=17, fontweight='bold',
            bbox=dict(facecolor='white', alpha=0.6, boxstyle='round,pad=0.25'))

# Annotate std
std = df[[A, T, C, G]].std()
min = df[[A, T, C, G]].min().min()
for i, val in enumerate(std):
    ax.text(i, min, f"std:{val:.1f}", ha='center', va='center',
            color='black', fontsize=17, fontweight='bold',
            bbox=dict(facecolor='gray', alpha=0.6, boxstyle='round,pad=0.25'))  '''


# Calculate mean and std
means = df[[A, T, C, G]].mean()
stds  = df[[A, T, C, G]].std()

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
plt.savefig('Chloro_N_usage_violon.svg', dpi=900, format='svg',transparent=True)
plt.show()