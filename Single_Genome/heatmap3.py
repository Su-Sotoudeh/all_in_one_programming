import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize

# 1. Data
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

# 2. DataFrame
df = pd.DataFrame(data, index=row_labels)
df.loc['by-'] = by_minus
df = df.loc[['CDS+', 'UTR+', 'Intron+', 'LINE+', 'SINE+', 'LTR+', 'DNA+', 'by-']]
df = df[['by+'] + list(reversed(['CDS-', 'UTR-', 'Intron-', 'LINE-', 'SINE-', 'LTR-', 'DNA-']))]

# 3. Custom colormap: light blue → soft red → light blue
colors = [
    (0.3, 0.4, 1.0),   # light blue (0.9, 1.1)
    (1.0, 0.3, 0.3),   # soft red (1.0)
    (0.3, 0.4, 1.0),
]
custom_cmap = LinearSegmentedColormap.from_list("soft_red_center", colors, N=256)
norm = Normalize(vmin=0.9, vmax=1.1)

# 4. Plot
fig, ax = plt.subplots(figsize=(5  , 4), dpi=200)
im = ax.imshow(df.values, cmap=custom_cmap, norm=norm)

# 5. Annotate
for i in range(df.shape[0]):
    for j in range(df.shape[1]):
        val = df.iloc[i, j]
        color = 'white' if abs(val - 1) > 0.05 else 'black'
        ax.text(j, i, f"{val:.3f}", ha='center', va='center', color=color, fontsize=6)

# 6. Axes and ticks
ax.set_xticks(np.arange(df.shape[1]))
ax.set_xticklabels(df.columns, rotation=45, ha='right', fontsize=8)
ax.set_yticks(np.arange(df.shape[0]))
ax.set_yticklabels(df.index, fontsize=8)
ax.tick_params(left=False, bottom=False)
for spine in ax.spines.values():
    spine.set_visible(False)

# 7. Colorbar
cbar = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
cbar.set_label('Ratio (normalized)', rotation=270, labelpad=10, fontsize=9)
cbar.ax.tick_params(labelsize=7.5)

# 8. Title
plt.title("Genomic Ratio Deviation (Red = Balanced, Blue = Divergent)", fontsize=9.5)
plt.tight_layout()
plt.show()
