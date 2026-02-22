import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize

# 1. Load your data
data = {
    'by+':      [1.177817526, 0.933379517, 0.911823548, 1.394453793, 1.384551574, 0.981038728, 1.015925344],
    'CDS-':     [1.007762939, 0.907269221, 0.909496373, 1.360091935, 1.330324631, 0.964416574, 0.974832608],
    'UTR-':     [1.102928358, 0.994909998, 0.922494834, 1.358158392, 1.329637765, 0.999105284, 1.034830132],
    'Intron-':  [1.102370063, 1.083445829, 0.996570659, 1.223310025, 1.185318103, 1.069608975, 1.088360534],
    'LINE-':    [0.735463109, 0.736006913, 0.822573181, 0.999730159, 0.916589378, 0.78134783, 0.753823965],
    'SINE-':    [0.753378554, 0.752748361, 0.84866839, 1.089917599, 0.999443974, 0.813892047, 0.77935533],
    'LTR-':     [1.035245377, 0.995150075, 0.933423754, 1.278656321, 1.227857189, 0.997784193, 1.015411881],
    'DNA-':     [1.028170397, 0.963582556, 0.919921953, 1.323240409, 1.280288978, 0.981996359, 0.999685544],
}
row_labels = ['CDS+', 'UTR+', 'Intron+', 'LINE+', 'SINE+', 'LTR+', 'DNA+']
by_minus = [1.0, 0.852290595, 1.065581406, 1.098979203, 0.715776433, 0.721623033, 1.015176126, 0.984078191]

# 2. Create DataFrame
df = pd.DataFrame(data, index=row_labels)
df.loc['by-'] = by_minus

# 3. Reorder rows and columns
df = df.loc[['CDS+', 'UTR+', 'Intron+', 'LINE+', 'SINE+', 'LTR+', 'DNA+', 'by-']]
reversed_cols = ['by+'] + list(reversed(['CDS-', 'UTR-', 'Intron-', 'LINE-', 'SINE-', 'LTR-', 'DNA-']))
df = df[reversed_cols]

# 4. Build custom colormap: Blue for far from 1, Red for 1
colors = [
    (0, 0, 1),       # blue at min
    (0.8, 0.8, 1),   # light blue
    (1, 0, 0),       # red at 1
    (0.8, 0.8, 1),   # light blue
    (0, 0, 1),       # blue at max
]
custom_cmap = LinearSegmentedColormap.from_list("custom_diverging", colors, N=256)

# 5. Normalize: center at 1, blue far from it
vmin, vmax = 0.90, 1.1  # min and max for color scale
norm = Normalize(vmin=vmin, vmax=vmax)

# 6. Plot
fig, ax = plt.subplots(figsize=(11, 7))
im = ax.imshow(df.values, cmap=custom_cmap, norm=norm)

# Annotate each cell
for i in range(df.shape[0]):
    for j in range(df.shape[1]):
        val = df.iloc[i, j]
        color = 'white' if abs(val - 1) > 0.02 else 'black'
        ax.text(j, i, f"{val:.3f}", ha='center', va='center', color=color, fontsize=9)

# Ticks and labels
ax.set_xticks(np.arange(df.shape[1]))
ax.set_xticklabels(df.columns, rotation=45, ha='right')
ax.set_yticks(np.arange(df.shape[0]))
ax.set_yticklabels(df.index)

# Cleanup
ax.tick_params(left=False, bottom=False)
for spine in ax.spines.values():
    spine.set_visible(False)

# Colorbar
cbar = fig.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label('Ratio value', rotation=270, labelpad=15)
cbar.ax.tick_params(labelsize=9)

# Title and layout
plt.title("Genomic Heatmap with Color Meaning", fontsize=13)
plt.tight_layout()
plt.show()
