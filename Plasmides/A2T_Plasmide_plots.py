

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# =========================================================
# LOAD DATA
# =========================================================

df = pd.read_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\Plasmides\Plasmides_refrence_Int_UTR_CDS_2026-02-05_23-22-23\Plasmides_2026-02-05_23-22-23_add_columns.csv')


# Optional sanity filter (same philosophy as before)
#df = df[(df['A2T_W'] > 0.8) & (df['A2T_W'] < 1.2)]

# =========================================================
# SINGLE CONDITION: P != 0 AND C != 0
# =========================================================

sub_df = df[
    (df['annot_P_A'] != 0) &
    (df['annot_C_A'] != 0)
    ]

# =========================================================
# FIGURE SETUP (2 rows × 1 column)
# =========================================================

fig, axes = plt.subplots(
    nrows=2,
    ncols=1,
    figsize=(8, 8),
    sharex=True
)

COLOR_P = '#008B8B'
COLOR_C = '#D55E00'
COLOR_A = 'green'
COLOR_W = '#333333'

size = 1.2
YMIN, YMAX = 0, 2

# =========================================================
# A2T ROW
# =========================================================

ax1 = axes[0]
dfa = sub_df.sort_values('A2T_P-A2T_C')
x = np.arange(len(dfa))

ax1.scatter(x, dfa['A2T_P'], color=COLOR_P, s=size)
ax1.scatter(x, dfa['A2T_C'], color=COLOR_C, s=size)
ax1.scatter(x, dfa['A2T_annot'], color=COLOR_A, s=size)
ax1.scatter(x, dfa['A2T_W'], color=COLOR_W, s=size)

ax1.axhline(1, linestyle='--', linewidth=1, alpha=0.7)
ax1.set_ylim(YMIN, YMAX)
ax1.set_xlim(0, len(x))
ax1.set_ylabel('A2T Ratio', fontsize=14, fontweight='bold')
ax1.set_title('Plasmids', fontsize=14, fontweight='bold')
ax1.tick_params(axis='x', labelbottom=False)

# ---- A2T STATS TABLE ----
stats_A2T = [
    ['A/T(+)',  f"{dfa['A2T_P'].mean():.3f}",     f"{dfa['A2T_P'].std():.3f}"],
    ['A/T(-)',  f"{dfa['A2T_C'].mean():.3f}",     f"{dfa['A2T_C'].std():.3f}"],
    ['A/T(+-)', f"{dfa['A2T_annot'].mean():.3f}", f"{dfa['A2T_annot'].std():.3f}"],
    ['Whole',   f"{dfa['A2T_W'].mean():.3f}",     f"{dfa['A2T_W'].std():.3f}"],
]

table = ax1.table(
    cellText=stats_A2T,
    colLabels=['', 'Mean', 'Std'],
    cellLoc='center',
    loc='bottom',
    bbox=[0, 0.0, 1, 0.28]
)
table.auto_set_font_size(False)
table.set_fontsize(9)

# =========================================================
# G2C ROW
# =========================================================

ax2 = axes[1]
dfg = sub_df.sort_values('G2C_P-G2C_C')
x = np.arange(len(dfg))

ax2.scatter(x, dfg['G2C_P'], color=COLOR_P, s=size)
ax2.scatter(x, dfg['G2C_C'], color=COLOR_C, s=size)
ax2.scatter(x, dfg['G2C_annot'], color=COLOR_A, s=size)
ax2.scatter(x, dfg['G2C_W'], color=COLOR_W, s=size)

ax2.axhline(1, linestyle='--', linewidth=1, alpha=0.7)
ax2.set_ylim(YMIN, YMAX)
ax2.set_xlim(0, len(x))
ax2.set_ylabel('G2C Ratio', fontsize=14, fontweight='bold')
ax2.tick_params(axis='x', labelsize=9)

# ---- G2C STATS TABLE ----
stats_G2C = [
    ['G/C(+)',  f"{dfg['G2C_P'].mean():.3f}",     f"{dfg['G2C_P'].std():.3f}"],
    ['G/C(-)',  f"{dfg['G2C_C'].mean():.3f}",     f"{dfg['G2C_C'].std():.3f}"],
    ['G/C(+-)', f"{dfg['G2C_annot'].mean():.3f}", f"{dfg['G2C_annot'].std():.3f}"],
    ['Whole',   f"{dfg['G2C_W'].mean():.3f}",     f"{dfg['G2C_W'].std():.3f}"],
]

table = ax2.table(
    cellText=stats_G2C,
    colLabels=['', 'Mean', 'Std'],
    cellLoc='center',
    loc='bottom',
    bbox=[0, 0.0, 1, 0.28]
)
table.auto_set_font_size(False)
table.set_fontsize(9)

# =========================================================
# LEGEND
# =========================================================

legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='+Annotation', markerfacecolor=COLOR_P, markersize=8),
    Line2D([0], [0], marker='o', color='w', label='-Annotation', markerfacecolor=COLOR_C, markersize=8),
    Line2D([0], [0], marker='o', color='w', label='+-Annotation', markerfacecolor=COLOR_A, markersize=8),
    Line2D([0], [0], marker='o', color='w', label='Whole', markerfacecolor=COLOR_W, markersize=8),
]

fig.legend(
    handles=legend_elements,
    loc='lower center',
    ncol=4,
    fontsize=12,
    bbox_to_anchor=(0.5, -0.02)
)

plt.tight_layout(rect=[0, 0.08, 1, 1])
plt.show()
