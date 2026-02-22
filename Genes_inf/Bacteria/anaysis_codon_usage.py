import matplotlib
matplotlib.use("Agg")   # non-GUI backend
import matplotlib.pyplot as plt

file_path = r"/home/genomics/Genomics/Coding/all_in_one_programming/Bacteria/codon_usage_df_6000.csv"
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
df = pd.read_csv(file_path)
df = df[df['GC%']>20]
print(df.head())
GC = df['GC%']
AAA = df['AAA']

plt.scatter(GC,df['ACA'],color='g',s= 0.5)
plt.scatter(GC,df['ACG'],color='r',s= 0.5)
plt.scatter(GC,df['ACT'],color='b',s= 0.5)
plt.scatter(GC,df['ACC'],color='y',s= 0.5)
plt.savefig("GC.png", dpi=600)



import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Load your dataframe (already created in your pipeline)
df = pd.read_csv("codon_usage_df.csv")

# Define codons (all columns except Organism and GC%)
codons = [c for c in df.columns if c not in ["Organism", "GC%"]]

# === Bin GC% into ranges ===
bins = np.arange(20, 80, 1)  # bins from 25% to 75% GC in steps of 5
df["GC_bin"] = pd.cut(df["GC%"], bins)

# === Compute average codon frequencies per bin ===
heatmap_data = df.groupby("GC_bin")[codons].mean().T  # codons as rows, GC bins as columns

# === Plot heatmap ===
plt.figure(figsize=(14, 18))
sns.heatmap(heatmap_data, cmap="viridis", cbar_kws={'label': 'Average frequency (%)'},vmax=10)
plt.xlabel("GC% bin")
plt.ylabel("Codon")
plt.title("Codon usage vs GC% bins")
plt.tight_layout()
plt.savefig("codon_usage_heatmap.png")
plt.close()






import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from Bio.Data import CodonTable

# Load data
df = pd.read_csv("codon_usage_df_6000.csv")

# GC% bins
bins = np.arange(25, 80, 5)
df["GC_bin"] = pd.cut(df["GC%"], bins)

# Genetic code
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
codon_to_aa = {c:a for c,a in standard_table.forward_table.items()}
for c in standard_table.stop_codons:
    codon_to_aa[c] = "Stop"

# Codon list
codons = [c for c in df.columns if c not in ["Organism", "GC%", "GC_bin"]]

# Build AA grouping
aa = [codon_to_aa.get(c, "Unknown") for c in codons]
aa_series = pd.Series(aa, index=codons, name="AA")

# Sort amino acids by GC-rich prefix score, then codons within AA
def gc_prefix_score(codon_list):
    return -sum(c[:2].count("G") + c[:2].count("C") for c in codon_list)

aa_to_codons = {}
for c in codons:
    aa_to_codons.setdefault(codon_to_aa.get(c, "Unknown"), []).append(c)

sorted_aa = sorted(aa_to_codons.keys(), key=lambda k: gc_prefix_score(aa_to_codons[k]))
codons_sorted = []
for k in sorted_aa:
    codons_sorted.extend(sorted(aa_to_codons[k]))

# Average per GC bin, rows=codons
heatmap_data = df.groupby("GC_bin")[codons_sorted].mean().T

# Row colors by AA
palette = sns.color_palette("tab20", n_colors=len(sorted_aa))
aa_color_map = {aa_name: palette[i] for i, aa_name in enumerate(sorted_aa)}
row_colors = heatmap_data.index.to_series().map(lambda c: aa_color_map[codon_to_aa.get(c, "Unknown")])

# Plot clustermap without clustering (we keep your order)
g = sns.clustermap(
    heatmap_data,
    row_cluster=False,
    col_cluster=False,
    cmap="coolwarm",
    vmin=0, vmax=10,
    row_colors=row_colors,
    cbar_kws={'label': 'Average frequency (%)'},
    figsize=(12, 18)
)

# Add legend for amino acids
from matplotlib.patches import Patch
handles = [Patch(facecolor=aa_color_map[a], label=a) for a in sorted_aa]
g.ax_heatmap.legend(handles=handles, title="Amino acid", loc='upper right', bbox_to_anchor=(1.25, 1))

g.ax_heatmap.set_xlabel("GC% bin")
g.ax_heatmap.set_ylabel("Codon (grouped by amino acid)")
g.fig.suptitle("Codon usage vs GC% bins (row-colored by amino acid)", y=1.02)
plt.savefig("codon_usage_clustermap_by_aa.png", bbox_inches="tight")
plt.close()


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from Bio.Data import CodonTable

df = pd.read_csv("codon_usage_df_6000.csv")
bins = np.arange(25, 80, 5)
df["GC_bin"] = pd.cut(df["GC%"], bins)

standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
codon_to_aa = {c:a for c,a in standard_table.forward_table.items()}
for c in standard_table.stop_codons:
    codon_to_aa[c] = "Stop"

codons = [c for c in df.columns if c not in ["Organism", "GC%", "GC_bin"]]

# Group and sort
aa_to_codons = {}
for c in codons:
    aa_to_codons.setdefault(codon_to_aa.get(c, "Unknown"), []).append(c)

def gc_prefix_score(codon_list):
    return -sum(c[:2].count("G") + c[:2].count("C") for c in codon_list)

sorted_aa = sorted(aa_to_codons.keys(), key=lambda k: gc_prefix_score(aa_to_codons[k]))
codons_sorted = []
aa_labels_per_row = []
for a in sorted_aa:
    for c in sorted(aa_to_codons[a]):
        codons_sorted.append(c)
        aa_labels_per_row.append(a)

heatmap_data = df.groupby("GC_bin")[codons_sorted].mean().T  # rows=codons_sorted

plt.figure(figsize=(12, 18))
ax = sns.heatmap(
    heatmap_data,
    cmap="coolwarm",
    vmin=0, vmax=10,
    cbar_kws={'label': 'Average frequency (%)'}
)

# Draw separators between amino acid groups
boundaries = []
offset = 0
for a in sorted_aa:
    count = len(sorted(aa_to_codons[a]))
    offset += count
    boundaries.append(offset)

# Lines at boundaries: use row centers with +/-0.5 offset
for b in boundaries[:-1]:
    ax.axhline(b, xmin=0, xmax=heatmap_data.shape[1], color="black", linewidth=1)
# Important: axhline uses x in data units only if transform is default; we need to span full width:
ax.set_xlim(0, heatmap_data.shape[1])

# Right-side secondary y-axis with amino acid labels per group (one label placed at group center)
ax2 = ax.twinx()
ax2.set_ylim(ax.get_ylim())  # align y-axis
group_centers = []
start = 0
for a in sorted_aa:
    count = len(sorted(aa_to_codons[a]))
    center = start + count / 2
    group_centers.append((center, a))
    start += count

ax2.set_yticks([c for c, _ in group_centers])
ax2.set_yticklabels([a for _, a in group_centers])
ax2.tick_params(axis='y', which='both', length=0)
ax2.set_ylabel("Amino acid groups")

ax.set_xlabel("GC% bin")
ax.set_ylabel("Codon")
ax.set_title("Codon usage vs GC% bins (grouped by amino acid with separators)")
plt.tight_layout()
plt.savefig("codon_usage_heatmap_grouped_with_separators.png")
plt.close()


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from Bio.Data import CodonTable

# === Load codon usage dataframe ===
df = pd.read_csv("codon_usage_df_6000.csv")

# === Define GC% bins ===
bins = np.arange(24, 80, 2)  # 25–75% GC in steps of 5
df["GC_bin"] = pd.cut(df["GC%"], bins)

# === Codon → Amino Acid mapping ===
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
codon_to_aa = {c: a for c, a in standard_table.forward_table.items()}
for c in standard_table.stop_codons:
    codon_to_aa[c] = "Stop"

# === Get codons from dataframe ===
codons = [c for c in df.columns if c not in ["Organism", "GC%", "GC_bin"]]

# === Compute average codon frequencies per GC bin ===
heatmap_data = df.groupby("GC_bin")[codons].mean().T

# === Build row colors for amino acid groups ===
aa_list = [codon_to_aa.get(c, "Unknown") for c in heatmap_data.index]
unique_aa = sorted(set(aa_list))
palette = sns.color_palette("tab20", n_colors=len(unique_aa))
aa_color_map = {aa: palette[i] for i, aa in enumerate(unique_aa)}
row_colors = pd.Series(aa_list, index=heatmap_data.index).map(aa_color_map)

# === Clustered heatmap with dendrograms ===
sns.set(style="white")
g = sns.clustermap(
    heatmap_data,
    cmap="viridis",
    vmin=0, vmax=10,
    row_cluster=True,   # cluster codons
    col_cluster=True,   # cluster GC bins
    row_colors=row_colors,
    cbar_kws={'label': 'Average frequency (%)'},
    figsize=(12, 18)
)

# === Add legend for amino acids ===
from matplotlib.patches import Patch
handles = [Patch(facecolor=aa_color_map[a], label=a) for a in unique_aa]
g.ax_heatmap.legend(handles=handles, title="Amino acid",
                    loc='upper right', bbox_to_anchor=(1.25, 1))

# === Labels ===
g.ax_heatmap.set_xlabel("GC% bin")
g.ax_heatmap.set_ylabel("Codons (clustered)")
g.fig.suptitle("Clustered heatmap of codon usage vs GC% bins", y=1.02)

# === Save ===
plt.savefig("codon_usage_clustermap.png", bbox_inches="tight")
plt.close()





