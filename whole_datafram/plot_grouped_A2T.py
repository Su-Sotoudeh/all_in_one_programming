import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Load data
df = pd.read_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\whole_datafram\Whole_datafram_2026.csv', low_memory=False)

# ==========================================
# 1. CORRECT GROUP NAMES (Scientific Names)
# ==========================================
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

# Apply the name changes to the DataFrame
df['Group'] = df['Group'].replace(name_mapping)

# ==========================================
# 2. PREFIX FILTERING (Your Logic)
# ==========================================
# Extract the first 2 or 3 characters to be safe
df['Prefix'] = df['Chr_ID'].str[:2]

# Define "Good" prefixes
good_prefixes = ['AC','AE','AL','AM','AP','AY','BA','BX','CM','CP','CR','FN','FR','GK','HE','HF','HG','L7','LN','LR','LT','NC','NZ','OU','OV','OW','OX','OY','OZ']

# Filter: Keep good prefixes
df_filtered = df[df['Prefix'].isin(good_prefixes)].copy()

print(f"Original shape: {df.shape}")
print(f"Filtered shape: {df_filtered.shape}")

# Update df to use the filtered version
df = df_filtered
df = df[(df['A2T_W']>0.9) & (df['A2T_W']<1.1)]
# ==========================================
# 3. PLOTTING LOOP (Rest of your code)
# ==========================================
groups = df['Group'].unique()
target_groups = groups[:8]

# Create figure
fig, axes = plt.subplots(nrows=2, ncols=8, figsize=(24, 8), sharey=True)
ax_flat = axes.flatten()

# Professional Colors
COLOR_P = '#008B8B'    # Teal
COLOR_C = '#D55E00'    # Vermilion
COLOR_A = '#0072B2'    # Blue
COLOR_A = 'green'
COLOR_W = '#333333'    # Dark Grey

size = 0.5

for i, group in enumerate(target_groups):
    dfg = df[df['Group'] == group].copy()

    # print(f"Group: {group}") # Optional debug
    # print(dfg['Chr_ID'].str[:3].unique())

    if dfg.empty:
        continue

    # --- ROW 1: A2T ---
    ax1 = ax_flat[i]
    dfg_sorted_A2T = dfg.sort_values('A2T_P-A2T_C')
    x = np.arange(len(dfg_sorted_A2T))
    sample_count = len(x)

    ax1.scatter(x, dfg_sorted_A2T['A2T_P'], color=COLOR_P, s=size)
    ax1.scatter(x, dfg_sorted_A2T['A2T_C'], color=COLOR_C, s=size)
    ax1.scatter(x, dfg_sorted_A2T['A2T_annot'], color=COLOR_A, s=size)
    ax1.scatter(x, dfg_sorted_A2T['A2T_W'], color=COLOR_W, s=size)

    ax1.set_ylim(0.2, 1.5)
    ax1.set_xlim(0, sample_count)

    # Titles: Group Name (Now uses Scientific Names)
    if i < 8:
        ax1.set_title(group, fontsize=12, fontweight='bold')
    else:
        ax1.set_title('')

    # Y-Axis: Ratio ONLY for the first column
    if i == 0:
        ax1.set_ylabel('A2T Ratio', fontsize=14, fontweight='bold')
        ax1.tick_params(axis='y', labelleft=True, labelsize=10)
    else:
        ax1.set_ylabel('')
        ax1.tick_params(axis='y', labelleft=False)

    ax1.tick_params(axis='x', labelbottom=False)

    # --- ADD A2T TABLE ---
    stats_A2T = [
        ['A/T(+)', f"{dfg['A2T_P'].mean():.3f}", f"{dfg['A2T_P'].std():.3f}"],
        ['A/T(-)', f"{dfg['A2T_C'].mean():.3f}", f"{dfg['A2T_C'].std():.3f}"],
        ['A/T(+-)', f"{dfg['A2T_annot'].mean():.3f}", f"{dfg['A2T_annot'].std():.3f}"],
        ['Whole', f"{dfg['A2T_W'].mean():.3f}", f"{dfg['A2T_W'].std():.3f}"]
    ]
    table = ax1.table(cellText=stats_A2T, colLabels=['', 'Mean', 'Std'],
                      cellLoc='center', loc='bottom', bbox=[0.0, 0.0, 1.0, 0.25])
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.2)

    # --- ROW 2: G2C ---
    ax2 = ax_flat[i + 8]
    dfg_sorted_G2C = dfg.sort_values('G2C_P-G2C_C')
    x = np.arange(len(dfg_sorted_G2C))
    sample_count = len(x)

    ax2.scatter(x, dfg_sorted_G2C['G2C_P'], color=COLOR_P, s=size)
    ax2.scatter(x, dfg_sorted_G2C['G2C_C'], color=COLOR_C, s=size)
    ax2.scatter(x, dfg_sorted_G2C['G2C_annot'], color=COLOR_A, s=size)
    ax2.scatter(x, dfg_sorted_G2C['G2C_W'], color=COLOR_W, s=size)

    ax2.set_ylim(0.2, 1.5)
    ax2.set_xlim(0, sample_count)

    # Y-Axis: Ratio ONLY for the first column
    if i == 0:
        ax2.set_ylabel('G2C Ratio', fontsize=14, fontweight='bold')
        ax2.tick_params(axis='y', labelleft=True, labelsize=10)
    else:
        ax2.set_ylabel('')
        ax2.tick_params(axis='y', labelleft=False)

    ax2.tick_params(axis='x', labelbottom=True, labelsize=9)

    # --- ADD G2C TABLE ---
    stats_G2C = [
        ['G/C(+)', f"{dfg['G2C_P'].mean():.3f}", f"{dfg['G2C_P'].std():.3f}"],
        ['G/C(-)', f"{dfg['G2C_C'].mean():.3f}", f"{dfg['G2C_C'].std():.3f}"],
        ['G/C(+-)', f"{dfg['G2C_annot'].mean():.3f}", f"{dfg['G2C_annot'].std():.3f}"],
        ['Whole', f"{dfg['G2C_W'].mean():.3f}", f"{dfg['G2C_W'].std():.3f}"]
    ]
    table = ax2.table(cellText=stats_G2C, colLabels=['', 'Mean', 'Std'],
                      cellLoc='center', loc='bottom', bbox=[0.0, 0.0, 1.0, 0.25])
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.2)

# Adjust layout
plt.tight_layout(rect=[0, 0.05, 1, 1])

# Add a single legend for the whole figure
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='+Annotation)', markerfacecolor=COLOR_P, markersize=8),
    Line2D([0], [0], marker='o', color='w', label='-Annotation', markerfacecolor=COLOR_C, markersize=8),
    Line2D([0], [0], marker='o', color='w', label='+-Annotation', markerfacecolor=COLOR_A, markersize=8),
    Line2D([0], [0], marker='o', color='w', label='Whole sequence', markerfacecolor=COLOR_W, markersize=8)
]
fig.legend(handles=legend_elements, loc='lower center', ncol=4, fontsize=12, bbox_to_anchor=(0.5, -0.02))

plt.show()

# ==========================================
# 4. STATISTICS TABLE (Prefix Analysis)
# ==========================================
# 1. Define the function to extract prefix
def get_prefix(id_str):
    s = str(id_str)[:3]
    if len(s) > 2 and s[2].isdigit():
        return s[:2]
    else:
        return s

# Apply the function
df['Prefix'] = df['Chr_ID'].apply(get_prefix)

# 2. Calculate Statistics
grouped = df.groupby(['Group', 'Prefix'])['Length(M)'].agg(['mean', 'count'])
grouped = grouped.rename(columns={'mean': 'Avg_Length(M)', 'count': 'Frequency'})
stats_df = grouped.reset_index()

# 3. Pivot the table for better display
pivot_table = stats_df.pivot(index='Group', columns='Prefix')
pivot_table['Avg_Length(M)'] = pivot_table['Avg_Length(M)'].round(2)

# Display the table
print("Statistics by Group and Prefix:")
print(pivot_table)

# 4. Save to CSV
pivot_table.to_csv('Prefix_Stats_Full.csv')
print("\nFull table saved to 'Prefix_Stats_Full.csv'")