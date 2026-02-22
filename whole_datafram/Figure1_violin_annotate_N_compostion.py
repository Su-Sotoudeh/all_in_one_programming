import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# ======================================
# Load master file
# ======================================
file_path = r'C:\Genomics\Chargaff\coding\all_in_one_programming\whole_datafram\Whole_datafram_2026_Eu_Pro_Virus_Plasmide.csv'
df_raw = pd.read_csv(file_path, low_memory=False)
# Convert all non-group columns to numeric safely
for col in df_raw.columns:
    if col != "Group":
        df_raw[col] = pd.to_numeric(df_raw[col], errors="coerce")

groups = df_raw['Group'].unique()
n_groups = len(groups)

# ======================================
# Biological Color Philosophy
# ======================================

COLOR_P = '#C62828'   # + annotation (Red)
COLOR_C = '#1565C0'   # − annotation (Blue)
COLOR_W = '#2E7D32'   # Combined ± (Green)
COLOR_A = '#000000'   # Chargaff rule (Black)
COLOR_U = '#7D7D7D'   # Unannotated (Gray)

row1_palette = [COLOR_P, COLOR_C, COLOR_W]

row2_palette = ['#FFCDD2','#EF9A9A','#E57373','#C62828']   # Red theme
row3_palette = ['#BBDEFB','#90CAF9','#64B5F6','#1565C0']   # Blue theme
row4_palette = ['#C8E6C9','#A5D6A7','#81C784','#2E7D32']   # Green theme
row5_palette = ['#C9A227', '#B8860B']

row6_palette = ['#E0E0E0','#BDBDBD','#9E9E9E','#616161']   # Gray theme
row7_palette = ['#FFD700', '#B8860B']  # Bright gold + dark gold

# ======================================
# Builders
# ======================================

def build_row1(df):
    P_anno = df['Int_P%'] + df['UTR_P%'] + df['CDS_P%'] + df['Share_P%']
    C_anno = df['Int_C%'] + df['UTR_C%'] + df['CDS_C%'] + df['Share_C%']
    Wh_anno = df['Wh_ano']

    return pd.DataFrame({
        'P_anno': P_anno,
        'C_anno': C_anno,
        'Wh_anno': Wh_anno
    })

def build_row2(df):
    return df[['annot_P_A','annot_P_T','annot_P_C','annot_P_G']]

def build_row3(df):
    return df[['annot_C_A','annot_C_T','annot_C_C','annot_C_G']]

def build_row4(df):
    P_size = df['Int_P%'] + df['UTR_P%'] + df['CDS_P%'] + df['Share_P%']
    C_size = df['Int_C%'] + df['UTR_C%'] + df['CDS_C%'] + df['Share_C%']

    P_portion = P_size/(P_size+C_size)
    C_portion = C_size/(P_size+C_size)

    result = pd.DataFrame()

    for base in ['A','T','C','G']:
        P_col = f'annot_P_{base}'
        C_col = f'annot_C_{base}'
        result[base] = ((df[P_col] * P_portion) + (df[C_col]*C_portion))

    return result

def build_row5(df):
    return df[['A2T_annot','G2C_annot']]

def build_row6(df):
    return df[['Unanot_P_A_ano','Unanot_P_T_ano',
               'Unanot_P_C_ano','Unanot_P_G_ano']]



def build_row7(df):
    return df[['Unanot_A2T_ano','Unanot_G2C_ano']]


# ======================================
# Scientific violin function
# ======================================

def plot_violin_with_stats(ax, data, palette, ylim):

    data = data.clip(lower=0)

    sns.violinplot(
        data=data,
        ax=ax,
        palette=palette,
        inner=None,
        cut=0,
        bw_adjust=0.7,
        scale='width',
        linewidth=1
    )

    sns.boxplot(
        data=data,
        ax=ax,
        width=0.25,
        boxprops={'facecolor':'none', 'edgecolor':'black'},
        showcaps=False,
        showfliers=False,
        whiskerprops={'linewidth':1}
    )

    means = data.mean()
    stds  = data.std()

    for i, col in enumerate(data.columns):

        mean_val = means[col]
        std_val  = stds[col]

        ax.scatter(i, mean_val,
                   color='black',
                   s=25,
                   zorder=5)

        ax.errorbar(i,
                    mean_val,
                    yerr=std_val,
                    fmt='none',
                    ecolor='black',
                    capsize=3,
                    linewidth=1)

        ax.text(
            i,
            ylim[1] * 0.97,
            f"{mean_val:.2f}±{std_val:.2f}",
            ha='center',
            va='top',
            fontsize=11,
            rotation=45,
            bbox=dict(
                boxstyle="round,pad=0.25",
                facecolor="white",
                edgecolor="none",
                alpha=0.6
            )
        )

    ax.set_ylim(ylim)
    ax.grid(alpha=0)

# ======================================
# Plot setup
# ======================================

sns.set(style="whitegrid", font_scale=0.9)

fig, axes = plt.subplots(
    7, n_groups,
    figsize=(3*n_groups, 21),
    sharey='row'
)

if n_groups == 1:
    axes = axes.reshape(7,1)

# ======================================
# Plot loop
# ======================================

for col, group in enumerate(groups):

    df_group = df_raw[df_raw['Group'] == group]

    # Row 1
    df1 = build_row1(df_group)
    ax = axes[0, col]
    plot_violin_with_stats(ax, df1, row1_palette, (0,100))
    ax.set_xticklabels(["P","C","Wh"])

    # Row 2
    df2 = build_row2(df_group)
    ax = axes[1, col]
    plot_violin_with_stats(ax, df2, row2_palette, (0,60))
    ax.set_xticklabels(["A","T","C","G"])

    # Row 3
    df3 = build_row3(df_group)
    ax = axes[2, col]
    plot_violin_with_stats(ax, df3, row3_palette, (0,60))
    ax.set_xticklabels(["A","T","C","G"])

    # Row 4
    df4 = build_row4(df_group)
    ax = axes[3, col]
    plot_violin_with_stats(ax, df4, row4_palette, (0,60))
    ax.set_xticklabels(["A","T","C","G"])

    # Row 5 (Chargaff rule)
    df5 = build_row5(df_group)
    ax = axes[4, col]
    plot_violin_with_stats(ax, df5, row5_palette, (0,2))
    ax.axhline(1, color='black', linestyle='--', linewidth=1)
    ax.set_xticklabels(["A/T","G/C"])

    # Row 6 (Unannotated)
    df6 = build_row6(df_group)
    ax = axes[5, col]
    plot_violin_with_stats(ax, df6, row6_palette, (0,60))
    ax.set_xticklabels(["A","T","C","G"])

    # Row 7 (Unannotated Chargaff)
    df7 = build_row7(df_group)
    ax = axes[6, col]
    plot_violin_with_stats(ax, df7, row7_palette, (0,2))
    ax.axhline(1, color='#B8860B', linestyle='--', linewidth=1.2)
    ax.set_xticklabels(["A/T","G/C"])





    axes[0, col].set_title(group, weight="bold")

# Row labels
axes[0,0].set_ylabel("Annotation Portion (%)")
axes[1,0].set_ylabel("+ Annotated (%)")
axes[2,0].set_ylabel("− Annotated (%)")
axes[3,0].set_ylabel("± Combined Annotation (%)")
axes[4,0].set_ylabel("Chargaff (Annotated)")
axes[5,0].set_ylabel("Unannotated Composition (%)")
axes[6,0].set_ylabel("Chargaff (Unannotated)")

import string

n_rows = axes.shape[0]
n_cols = axes.shape[1]

for i in range(n_rows):
    for j in range(n_cols):

        ax = axes[i, j]

        row_number = i + 1
        col_letter = string.ascii_uppercase[j]
        panel_label = f"{row_number}{col_letter}"

        ax.text(
            0.02, 0.02,              # bottom-left inside
            panel_label,
            transform=ax.transAxes,
            fontsize=13,
            fontweight='bold',
            ha='left',
            va='bottom'
        )






plt.tight_layout()
plt.savefig("Full_Composition_Matrix_6Rows_ColorLogic.svg",
            format="svg", dpi=900, transparent=True)

plt.show()
