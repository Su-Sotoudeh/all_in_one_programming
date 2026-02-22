import os
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# === File paths ===
Bacteria_file_path = r'C:\Genomics\results\Genes_chardaff_parameters\Bacterial_2025-09-23_20-16-24.csv'
Fungi_file_path = r'C:\Genomics\results\Genes_chardaff_parameters\Fungisresults_2025-09-17_00-53-54.csv'
Arthropoda_file_path = r'C:\Genomics\results\Genes_chardaff_parameters\Arthropoda_500m.csv'
Chordata_file_path = r'C:\Genomics\results\Genes_chardaff_parameters\Chordata.csv'
Viride_file_path = r'C:\Genomics\results\Genes_chardaff_parameters\Viride_on_saveing_time_resuts.csv'
Mito_file_path = r'C:\Genomics\results\Genes_chardaff_parameters\Mito_2025-11-18_19-16-56.csv'
# === Load CSVs ===
Bacteria_df = pd.read_csv(Bacteria_file_path, index_col=0)
Fungi_df = pd.read_csv(Fungi_file_path, index_col=0)
Arthropoda_df = pd.read_csv(Arthropoda_file_path, index_col=0)
Chordata_df = pd.read_csv(Chordata_file_path, index_col=0)
Viride_df = pd.read_csv(Viride_file_path, index_col=0)
Mito_df = pd.read_csv(Mito_file_path, index_col=0)
# === Column names ===
Wh_ano = 'Wh_ano'
A = 'CDS_P_A_ano'
T = 'CDS_P_T_ano'
C = 'CDS_P_C_ano'
G = 'CDS_P_G_ano'





df = Bacteria_df

# === STEP 4: Show statistics ===
print(f"\n--- {'Bacteria'} ---")
print(df.describe().to_string())


# === STEP 5: Plot ===
sns.set(style="whitegrid", font_scale=1.2)
fig, axes = plt.subplots(1, 4, figsize=(12, 6))
fig.suptitle(f"Nucleotide Composition Distributions - {'Bacteria'}", fontsize=14, weight="bold")

# ---------------- VIOLIN 0: A, T, C, G (with box inside) ----------------




ax = axes[0]
sns.violinplot(data=df[[Wh_ano]],
               ax=ax,
               palette=sns.light_palette("turquoise", n_colors=1),
               inner=None, linewidth=1.2, alpha=0.9)
sns.boxplot(data=df[[Wh_ano]],
            ax=ax,
            width=0.2, boxprops={'facecolor':'none', 'edgecolor':'black', 'linewidth':1.2},
            showcaps=False, showfliers=False, whiskerprops={'linewidth':0})
ax.set_title("Nucleotide Composition (Violin + Box)")
ax.set_ylabel("Percentage (%)")

means = df[[Wh_ano]].mean()
min =   df[[Wh_ano]].min()
for i, val in enumerate(means):
    ax.text(i-0.15, val, f"{val:.1f}", ha='center', va='center',
            color='black', fontsize=9, fontweight='bold',
            bbox=dict(facecolor='white', alpha=0.6, boxstyle='round,pad=0.25'))


sum_rations = df['Int_P%'] + df['UTR_P%'] + df['CDS_P%']
P2C = ((df['Int_P%']/sum_rations) * df['Int_P2C']) + ((df['UTR_P%']/sum_rations) * df['UTR_P2C']) + ((df['CDS_P%']/sum_rations) * df['CDS_P2C'])
ax = axes[1]
this_df = pd.DataFrame({'P2C': P2C})
sns.violinplot(data=this_df[['P2C']],
               ax=ax,
               palette=sns.light_palette("turquoise", n_colors=1),
               inner=None, linewidth=1.2, alpha=0.9)
sns.boxplot(data=this_df[['P2C']],
            ax=ax,
            width=0.2, boxprops={'facecolor':'none', 'edgecolor':'black', 'linewidth':1.2},
            showcaps=False, showfliers=False, whiskerprops={'linewidth':0})
ax.set_title("Nucleotide Composition (Violin + Box)")
ax.set_ylabel("Percentage (%)")

means = this_df[['P2C']].mean()
min =   this_df[['P2C']].min()
for i, val in enumerate(means):
    ax.text(i-0.15, val, f"{val:.1f}", ha='center', va='center',
            color='black', fontsize=9, fontweight='bold',
            bbox=dict(facecolor='white', alpha=0.6, boxstyle='round,pad=0.25'))



# ---------------- VIOLIN 1: A, T, C, G (with box inside) ----------------

Ac = 'CDS_C_A_ano'
Tc = 'CDS_C_T_ano'
Cc = 'CDS_C_C_ano'
Gc = 'CDS_C_G_ano'

ax = axes[2]
sns.violinplot(data=df[[A,T,C,G,Ac,Tc,Cc,Gc]],
               ax=ax,
               palette=sns.light_palette("turquoise", n_colors=8),
               inner=None, linewidth=1.2, alpha=0.9)
sns.boxplot(data=df[[A,T,C,G,Ac,Tc,Cc,Gc]],
            ax=ax,
            width=0.2, boxprops={'facecolor':'none', 'edgecolor':'black', 'linewidth':1.2},
            showcaps=False, showfliers=False, whiskerprops={'linewidth':0})
ax.set_title("Nucleotide Composition (Violin + Box)")
ax.set_ylabel("Percentage (%)")

means = df[[A,T,C,G,Ac,Tc,Cc,Gc]].mean()
min =   df[[A,T,C,G,Ac,Tc,Cc,Gc]].min()
for i, val in enumerate(means):
    ax.text(i-0.15, val, f"{val:.1f}", ha='center', va='center',
            color='black', fontsize=9, fontweight='bold',
            bbox=dict(facecolor='white', alpha=0.6, boxstyle='round,pad=0.25'))

# ---------------- VIOLIN 2: A/T, G/C (with box inside) ----------------

sum_rations = df['Int_P%'] + df['UTR_P%'] + df['CDS_P%']
A2T_P = ((df['Int_P%']/sum_rations) * df['Int_P_A2T']) + ((df['UTR_P%']/sum_rations) * df['UTR_P_A2T']) + ((df['CDS_P%']/sum_rations) * df['CDS_P_A2T'])
sum_rations = df['Int_P%'] + df['UTR_P%'] + df['CDS_P%']
G2C_P = ((df['Int_P%']/sum_rations) * df['Int_P_G2C']) + ((df['UTR_P%']/sum_rations) * df['UTR_P_G2C']) + ((df['CDS_P%']/sum_rations) * df['CDS_P_G2C'])

sum_rations = df['Int_C%'] + df['UTR_C%'] + df['CDS_C%']
A2T_C = ((df['Int_C%']/sum_rations) * df['Int_C_A2T']) + ((df['UTR_C%']/sum_rations) * df['UTR_C_A2T']) + ((df['CDS_C%']/sum_rations) * df['CDS_C_A2T'])
sum_rations = df['Int_C%'] + df['UTR_C%'] + df['CDS_C%']
G2C_C = ((df['Int_C%']/sum_rations) * df['Int_C_G2C']) + ((df['UTR_C%']/sum_rations) * df['UTR_C_G2C']) + ((df['CDS_C%']/sum_rations) * df['CDS_C_G2C'])



ratio_df = pd.DataFrame({'A2T_P': A2T_P,'G2C_P': G2C_P,'A2T_C': A2T_C,'G2C_C': G2C_C})

print(ratio_df.head())
ax = axes[3]
sns.violinplot(data=ratio_df[['A2T_P','G2C_P','A2T_C','G2C_C']],
               ax=ax,
               palette=sns.light_palette("teal", n_colors=4),
               inner=None, linewidth=1.2, alpha=0.9)
sns.boxplot(data=ratio_df[['A2T_P','G2C_P','A2T_C','G2C_C']],
            ax=ax,
            width=0.2, boxprops={'facecolor':'none', 'edgecolor':'black', 'linewidth':1.2},
            showcaps=False, showfliers=False, whiskerprops={'linewidth':0})
ax.set_title("A/T and G/C Ratios (Violin + Box)")
ax.set_ylabel("Ratio")

means = ratio_df[['A2T_P','G2C_P','A2T_C','G2C_C']].mean()
for i, val in enumerate(means):
    ax.text(i, val+0.1, f"{val:.2f}", ha='center', va='center',
            color='black', fontsize=9, fontweight='bold',
            bbox=dict(facecolor='white', alpha=0.6, boxstyle='round,pad=0.25'))




plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.show()


