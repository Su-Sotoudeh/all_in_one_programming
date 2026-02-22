

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

COLOR_P = '#C62828'# '#008B8B'    # Teal
COLOR_C = '#1565C0'#'#D55E00'    # Vermilion
COLOR_A = '#000000'#'green'      # Annot
COLOR_W = '#2E7D32'#'#333333'    # Whole
size = 0.5

YLIM_ROW1 = (0.4, 1.35)
YLIM_ROW2 = (0.3, 1.5)
YLIM_ROW3 = (0.0, 2.0)
YLIM_ROW4 = (0.0, 2.0)

# =========================================================
# LOAD DATASETS
# =========================================================

df1 = pd.read_csv(
    r'C:\Genomics\Chargaff\coding\all_in_one_programming\whole_datafram\Whole_datafram_2026.csv',
    low_memory=False
)

name_mapping = {
    'Chordata': 'Chordata',
    'Viride': 'Viridiplantae',
    'Insect': 'Arthropoda',
    'Fungi': 'Fungi',
    'Archia': 'Archaea',
    'Bacteria': 'Bacteria',
    'Chloro_plast': 'Chloroplast',
    'Mito': 'Mitochondria'
}
df1['Group'] = df1['Group'].replace(name_mapping)
df1['Prefix'] = df1['Chr_ID'].str[:2]

good_prefixes = [
    'AC','AE','AL','AM','AP','AY','BA','BX','CM','CP','CR','FN','FR','GK',
    'HE','HF','HG','L7','LN','LR','LT','NC','NZ','OU','OV','OW','OX','OY','OZ'
]
df1 = df1[df1['Prefix'].isin(good_prefixes)]
df1 = df1[(df1['A2T_W'] > 0.9) & (df1['A2T_W'] < 1.1)]

groups = df1['Group'].unique()
target_groups = groups[:8]

# ---------- VIRUSES ----------
v_main = pd.read_csv(
    r'C:\Genomics\Chargaff\coding\all_in_one_programming\Viruses\Viruses_refrence_Int_UTR_CDS_2026-01-26_16-51-15\Viruses_2026-01-26_16-51-15.csv',
    low_memory=False
)
v_ref = pd.read_csv(
    r"C:\Genomics\Genomes\unzip_viruse_all_complete\unzip_viruse_all_complete\Viruse_file_information\Refrence_viruses_list.csv"
)[['Assembly', 'Molecule_type']]

df2 = v_main.merge(v_ref, left_on='Name', right_on='Assembly', how='inner')
df2['Genome'] = df2['Molecule_type']

conditions = {
    '+':   (df2['annot_P_A'] != 0) & (df2['annot_C_A'] == 0),
    '-':   (df2['annot_P_A'] == 0) & (df2['annot_C_A'] != 0),
    '+/-': (df2['annot_P_A'] != 0) & (df2['annot_C_A'] != 0),
}
virus_layout = [
    ('RNA', '+'),
    ('RNA', '-'),
    ('RNA', '+/-'),
    ('DNA', '+'),
    ('DNA', '-'),
    ('DNA', '+/-'),
]

# ---------- PLASMIDS ----------
df3 = pd.read_csv(
    r'C:\Genomics\Chargaff\coding\all_in_one_programming\Plasmides\Plasmides_refrence_Int_UTR_CDS_2026-02-05_23-22-23\Plasmides_2026-02-05_23-22-23_add_columns.csv'
)
df3 = df3[(df3['annot_P_A'] != 0) & (df3['annot_C_A'] != 0)]

# ---------- ARTIFICIAL SIMULATION ----------
results_df = pd.read_csv(
    r'C:\Genomics\Chargaff\coding\all_in_one_programming\Artificial_Chromosome\simulation_results.csv'
)

# =========================================================
# CREATE MASTER FIGURE
# =========================================================
fig, axes = plt.subplots(4, 8, figsize=(38,25))
plt.subplots_adjust(hspace=0.55)

# =========================================================
# 1–2: EUKARYOTE + PROKARYOTE (Sheet 1 style)
# =========================================================
for i, group in enumerate(target_groups):
    dfg = df1[df1['Group'] == group]
    if dfg.empty:
        continue

    # ---------- A2T (Row 1) ----------
    ax = axes[0, i]
    dfa = dfg.sort_values('A2T_P-A2T_C')
    x = np.arange(len(dfa))

    ax.scatter(x, dfa['A2T_P'], color=COLOR_P, s=size)
    ax.scatter(x, dfa['A2T_C'], color=COLOR_C, s=size)
    ax.scatter(x, dfa['A2T_annot'], color=COLOR_A, s=size)
    ax.scatter(x, dfa['A2T_W'], color=COLOR_W, s=size)

    ax.axhline(1, linestyle='--', linewidth=1, alpha=0.7)
    ax.set_ylim(*YLIM_ROW1)
    ax.set_xlim(0, len(x))
    ax.set_title(group, fontsize=12, fontweight='bold')

    if i == 0:
        ax.set_ylabel("A2T Ratio", fontsize=14, fontweight='bold')
    ax.tick_params(axis='x', labelbottom=False)

    stats_A2T = [
        ['A/T(+)',  f"{dfa['A2T_P'].mean():.3f}", f"{dfa['A2T_P'].std():.3f}"],
        ['A/T(-)',  f"{dfa['A2T_C'].mean():.3f}", f"{dfa['A2T_C'].std():.3f}"],
        ['A/T(+-)', f"{dfa['A2T_annot'].mean():.3f}", f"{dfa['A2T_annot'].std():.3f}"],
        ['Whole',   f"{dfa['A2T_W'].mean():.3f}", f"{dfa['A2T_W'].std():.3f}"],
    ]
    table = ax.table(
        cellText=stats_A2T,
        colLabels=['Sequence', 'Mean', 'Std'],
        cellLoc='center',
        loc='bottom',
        bbox=[0.0, 0.0, 1.0, 0.25]
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(0.7, 1.2)
    #for cell in table.get_celld().values():
    #    cell.get_text().set_fontweight('bold')


    # ---------- G2C (Row 2) ----------
    ax = axes[1, i]
    dfg2 = dfg.sort_values('G2C_P-G2C_C')
    x = np.arange(len(dfg2))

    ax.scatter(x, dfg2['G2C_P'], color=COLOR_P, s=size)
    ax.scatter(x, dfg2['G2C_C'], color=COLOR_C, s=size)
    ax.scatter(x, dfg2['G2C_annot'], color=COLOR_A, s=size)
    ax.scatter(x, dfg2['G2C_W'], color=COLOR_W, s=size)

    ax.axhline(1, linestyle='--', linewidth=1, alpha=0.7)
    ax.set_ylim(*YLIM_ROW2)
    ax.set_xlim(0, len(x))

    if i == 0:
        ax.set_ylabel("G2C Ratio", fontsize=14, fontweight='bold')

    stats_G2C = [
        ['G/C(+)',  f"{dfg2['G2C_P'].mean():.3f}", f"{dfg2['G2C_P'].std():.3f}"],
        ['G/C(-)',  f"{dfg2['G2C_C'].mean():.3f}", f"{dfg2['G2C_C'].std():.3f}"],
        ['G/C(+-)', f"{dfg2['G2C_annot'].mean():.3f}", f"{dfg2['G2C_annot'].std():.3f}"],
        ['Whole',   f"{dfg2['G2C_W'].mean():.3f}", f"{dfg2['G2C_W'].std():.3f}"],
    ]
    table = ax.table(
        cellText=stats_G2C,
        colLabels=['Sequence', 'Mean', 'Std'],
        cellLoc='center',
        loc='bottom',
        bbox=[0.0, 0.0, 1.0, 0.25]
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.2)

# =========================================================
# 3–4: VIRUSES (first-sheet style)
# =========================================================
for col, (genome, cond) in enumerate(virus_layout):
    sub = df2[conditions[cond] & df2['Genome'].str.contains(genome, na=False)]
    if sub.empty:
        continue

    # ---------- A2T (Row 3) ----------
    ax = axes[2, col]
    dfa = sub.sort_values('A2T_P-A2T_C')
    x = np.arange(len(dfa))

    ax.scatter(x, dfa['A2T_P'], color=COLOR_P, s=size)
    ax.scatter(x, dfa['A2T_C'], color=COLOR_C, s=size)
    ax.scatter(x, dfa['A2T_annot'], color=COLOR_A, s=size)
    ax.scatter(x, dfa['A2T_W'], color=COLOR_W, s=size)

    ax.axhline(1, linestyle='--', linewidth=1, alpha=0.7)
    ax.set_ylim(*YLIM_ROW3)
    ax.set_xlim(0, len(x))
    ax.set_title(f"{cond}{genome} Viruses", fontsize=12, fontweight='bold')

    if col == 0:
        ax.set_ylabel("A2T Ratio", fontsize=14, fontweight='bold')
    ax.tick_params(axis='x', labelbottom=False)

    stats_A2T = [
        ['A/T(+)',  f"{dfa['A2T_P'].mean():.3f}", f"{dfa['A2T_P'].std():.3f}"],
        ['A/T(-)',  f"{dfa['A2T_C'].mean():.3f}", f"{dfa['A2T_C'].std():.3f}"],
        ['A/T(+-)', f"{dfa['A2T_annot'].mean():.3f}", f"{dfa['A2T_annot'].std():.3f}"],
        ['Whole',   f"{dfa['A2T_W'].mean():.3f}", f"{dfa['A2T_W'].std():.3f}"],
    ]
    table = ax.table(
        cellText=stats_A2T,
        colLabels=['Sequence', 'Mean', 'Std'],
        cellLoc='center',
        loc='bottom',
        bbox=[0.0, 0.0, 1.0, 0.25]
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.2)



    # ---------- G2C (Row 4) ----------
    ax = axes[3, col]
    dfg = sub.sort_values('G2C_P-G2C_C')
    x = np.arange(len(dfg))

    ax.scatter(x, dfg['G2C_P'], color=COLOR_P, s=size)
    ax.scatter(x, dfg['G2C_C'], color=COLOR_C, s=size)
    ax.scatter(x, dfg['G2C_annot'], color=COLOR_A, s=size)
    ax.scatter(x, dfg['G2C_W'], color=COLOR_W, s=size)

    ax.axhline(1, linestyle='--', linewidth=1, alpha=0.7)
    ax.set_ylim(*YLIM_ROW4)
    ax.set_xlim(0, len(x))

    if col == 0:
        ax.set_ylabel("G2C Ratio", fontsize=14, fontweight='bold')

    stats_G2C = [
        ['G/C(+)',  f"{dfg['G2C_P'].mean():.3f}", f"{dfg['G2C_P'].std():.3f}"],
        ['G/C(-)',  f"{dfg['G2C_C'].mean():.3f}", f"{dfg['G2C_C'].std():.3f}"],
        ['G/C(+-)', f"{dfg['G2C_annot'].mean():.3f}", f"{dfg['G2C_annot'].std():.3f}"],
        ['Whole',   f"{dfg['G2C_W'].mean():.3f}", f"{dfg['G2C_W'].std():.3f}"],
    ]
    table = ax.table(
        cellText=stats_G2C,
        colLabels=['Sequence', 'Mean', 'Std'],
        cellLoc='center',
        loc='bottom',
        bbox=[0.0, 0.0, 1.0, 0.25]
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.2)

# =========================================================
# 3–4: PLASMIDS (Column 7)
# =========================================================

# ---------- A2T ----------
ax = axes[2, 6]
dfa = df3.sort_values('A2T_P-A2T_C')
x = np.arange(len(dfa))

ax.scatter(x, dfa['A2T_P'], color=COLOR_P, s=size)
ax.scatter(x, dfa['A2T_C'], color=COLOR_C, s=size)
ax.scatter(x, dfa['A2T_annot'], color=COLOR_A, s=size)
ax.scatter(x, dfa['A2T_W'], color=COLOR_W, s=size)

ax.axhline(1, linestyle='--', linewidth=1, alpha=0.7)
ax.set_ylim(*YLIM_ROW3)
ax.set_xlim(0, len(x))
ax.set_title("Plasmids", fontsize=12, fontweight='bold')
ax.set_ylabel("A2T Ratio", fontsize=14, fontweight='bold')
ax.tick_params(axis='x', labelbottom=False)

stats_A2T = [
    ['A/T(+)',  f"{dfa['A2T_P'].mean():.3f}", f"{dfa['A2T_P'].std():.3f}"],
    ['A/T(-)',  f"{dfa['A2T_C'].mean():.3f}", f"{dfa['A2T_C'].std():.3f}"],
    ['A/T(+-)', f"{dfa['A2T_annot'].mean():.3f}", f"{dfa['A2T_annot'].std():.3f}"],
    ['Whole',   f"{dfa['A2T_W'].mean():.3f}", f"{dfa['A2T_W'].std():.3f}"],
]
table = ax.table(
    cellText=stats_A2T,
    colLabels=['Sequence', 'Mean', 'Std'],
    cellLoc='center',
    loc='bottom',
    bbox=[0.0, 0.0, 1.0, 0.25]
)
table.auto_set_font_size(False)
table.set_fontsize(8)
table.scale(1, 1.2)

# ---------- G2C ----------
ax = axes[3, 6]
dfg = df3.sort_values('G2C_P-G2C_C')
x = np.arange(len(dfg))

ax.scatter(x, dfg['G2C_P'], color=COLOR_P, s=size)
ax.scatter(x, dfg['G2C_C'], color=COLOR_C, s=size)
ax.scatter(x, dfg['G2C_annot'], color=COLOR_A, s=size)
ax.scatter(x, dfg['G2C_W'], color=COLOR_W, s=size)

ax.axhline(1, linestyle='--', linewidth=1, alpha=0.7)
ax.set_ylim(*YLIM_ROW4)
ax.set_xlim(0, len(x))
ax.set_ylabel("G2C Ratio", fontsize=14, fontweight='bold')

stats_G2C = [
    ['G/C(+)',  f"{dfg['G2C_P'].mean():.3f}", f"{dfg['G2C_P'].std():.3f}"],
    ['G/C(-)',  f"{dfg['G2C_C'].mean():.3f}", f"{dfg['G2C_C'].std():.3f}"],
    ['G/C(+-)', f"{dfg['G2C_annot'].mean():.3f}", f"{dfg['G2C_annot'].std():.3f}"],
    ['Whole',   f"{dfg['G2C_W'].mean():.3f}", f"{dfg['G2C_W'].std():.3f}"],
]
table = ax.table(
    cellText=stats_G2C,
    colLabels=['Sequence', 'Mean', 'Std'],
    cellLoc='center',
    loc='bottom',
    bbox=[0.0, 0.0, 1.0, 0.25]
)
table.auto_set_font_size(False)
table.set_fontsize(8)
table.scale(1, 1.2)

# =========================================================
# EMPTY CELLS (Row 3–4, Column 8) – turn off axes
# =========================================================



x_sim = np.arange(len(results_df))

# ---------- A2T Simulation (Row 3, Column 8) ----------
ax = axes[2, 7]

ax.scatter(x_sim, results_df['A2T_P'], color=COLOR_P, s=size)
ax.scatter(x_sim, results_df['A2T_C'], color=COLOR_C, s=size)
ax.scatter(x_sim, results_df['A2T_anno'], color=COLOR_A, s=size)

ax.axhline(1, linestyle='--', linewidth=1, alpha=0.7)
ax.set_ylim(*YLIM_ROW3)
ax.set_xlim(0, len(x_sim))
ax.set_title("Artificial Simulation", fontsize=12, fontweight='bold')
ax.set_ylabel("A2T Ratio", fontsize=14, fontweight='bold')
ax.tick_params(axis='x', labelbottom=False)

stats_A2T_sim = [
    ['A/T(+)',  f"{results_df['A2T_P'].mean():.3f}", f"{results_df['A2T_P'].std():.3f}"],
    ['A/T(-)',  f"{results_df['A2T_C'].mean():.3f}", f"{results_df['A2T_C'].std():.3f}"],
    ['A/T(+-)', f"{results_df['A2T_anno'].mean():.3f}", f"{results_df['A2T_anno'].std():.3f}"],
]

table = ax.table(
    cellText=stats_A2T_sim,
    colLabels=['Sequence', 'Mean', 'Std'],
    cellLoc='center',
    loc='bottom',
    bbox=[0.0, 0.0, 1.0, 0.25]
)
table.auto_set_font_size(False)
table.set_fontsize(8)
table.scale(1, 1.2)

# ---------- G2C Simulation (Row 4, Column 8) ----------
ax = axes[3, 7]

ax.scatter(x_sim, results_df['G2C_P'], color=COLOR_P, s=size)
ax.scatter(x_sim, results_df['G2C_C'], color=COLOR_C, s=size)
ax.scatter(x_sim, results_df['G2C_anno'], color=COLOR_A, s=size)

ax.axhline(1, linestyle='--', linewidth=1, alpha=0.7)
ax.set_ylim(*YLIM_ROW4)
ax.set_xlim(0, len(x_sim))
ax.set_ylabel("G2C Ratio", fontsize=14, fontweight='bold')

stats_G2C_sim = [
    ['G/C(+)',  f"{results_df['G2C_P'].mean():.3f}", f"{results_df['G2C_P'].std():.3f}"],
    ['G/C(-)',  f"{results_df['G2C_C'].mean():.3f}", f"{results_df['G2C_C'].std():.3f}"],
    ['G/C(+-)', f"{results_df['G2C_anno'].mean():.3f}", f"{results_df['G2C_anno'].std():.3f}"],
]

table = ax.table(
    cellText=stats_G2C_sim,
    colLabels=['Sequence', 'Mean', 'Std'],
    cellLoc='center',
    loc='bottom',
    bbox=[0.0, 0.0, 1.0, 0.25]
)
table.auto_set_font_size(False)
table.set_fontsize(8)
table.scale(1, 1.2)

# =========================================================
# GLOBAL LEGEND (UNCHANGED)
# =========================================================

legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='+Annotation', markerfacecolor=COLOR_P, markersize=8),
    Line2D([0], [0], marker='o', color='w', label='-Annotation', markerfacecolor=COLOR_C, markersize=8),
    Line2D([0], [0], marker='o', color='w', label='+-Annotation', markerfacecolor=COLOR_A, markersize=8),
    Line2D([0], [0], marker='o', color='w', label='Whole sequence', markerfacecolor=COLOR_W, markersize=8),
]

fig.legend(
    handles=legend_elements,
    loc='lower center',
    ncol=4,
    fontsize=16,
    bbox_to_anchor=(0.5, 0.01)
)

plt.tight_layout(rect=[0, 0.08, 1, 1])
plt.show()
