import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def plot_signature_violin(
        df,
        cell_type_col="cell_type",
        value_col="signature_expr",
        order=None,
        title="scRNA-seq huBCC",
        subtitle="Spatial Tumor DEG signature (57 genes)",
        ylabel="Expression level (in normalized counts +1)",
        palette="viridis",
        log1p=False,
        show_points=False,
        figsize=(12, 5),
):
    data = df[[cell_type_col, value_col]].dropna().copy()

    if log1p:
        data["_plot_val"] = np.log1p(data[value_col])
        ylab = f"{ylabel} (log1p)"
    else:
        data["_plot_val"] = data[value_col]
        ylab = ylabel

    if order is None:
        order = sorted(data[cell_type_col].unique())

    sns.set_theme(context="talk", style="white")
    fig, ax = plt.subplots(figsize=figsize)

    # ✅ Updated for Seaborn 0.14+ compatibility
    sns.violinplot(
        data=data,
        x=cell_type_col,
        y="_plot_val",
        hue=cell_type_col,
        order=order,
        inner="quartile",
        cut=0,
        density_norm="width",  # replaces scale="width"
        linewidth=1,
        palette=sns.color_palette(palette, n_colors=len(order)),
        legend=False,
        ax=ax
    )

    if show_points:
        sns.stripplot(
            data=data,
            x=cell_type_col,
            y="_plot_val",
            order=order,
            color="black",
            size=2,
            jitter=0.15,
            alpha=0.35,
            ax=ax
        )

    ax.set_title(subtitle, fontsize=14, pad=6)
    fig.suptitle(title, fontsize=18, y=0.98)
    ax.set_xlabel("")
    ax.set_ylabel(ylab)

    sns.despine(offset=5, trim=True)
    ax.grid(False)
    ax.tick_params(axis="x", rotation=20)
    plt.tight_layout()
    plt.show()


# ----------------------------
# Example usage with synthetic data
# ----------------------------

np.random.seed(42)
cell_types = [
    "Epithelial cells", "T cells", "Macrophages", "Myeloid cells",
    "Endothelial cells", "Fibroblasts", "Cycling cells", "Pericytes",
    "Mast cells", "B cells"
]

rows = []
for ct in cell_types:
    n = np.random.randint(120, 250)
    base = np.random.uniform(0.5, 2.0)
    vals = np.random.gamma(shape=2.0, scale=base, size=n) + 1.0
    rows.extend([(ct, v) for v in vals])

df_example = pd.DataFrame(rows, columns=["cell_type", "signature_expr"])

order = [
    "Epithelial cells", "T cells", "Macrophages", "Myeloid cells",
    "Endothelial cells", "Fibroblasts", "Cycling cells", "Pericytes",
    "Mast cells", "B cells"
]

plot_signature_violin(
    df_example,
    cell_type_col="cell_type",
    value_col="signature_expr",
    order=order,
    title="scRNA-seq huBCC",
    subtitle="Spatial Tumor DEG signature (57 genes)",
    ylabel="Expression level (in normalized counts +1)",
    palette="viridis",
    log1p=False,
    show_points=False,
    figsize=(12, 5),
)


import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Synthetic example data
np.random.seed(42)
cell_types = ["Epithelial", "T cells", "Macrophages", "Fibroblasts"]
rows = []
for ct in cell_types:
    vals = np.random.gamma(shape=2.0, scale=1.5, size=120) + 1.0
    rows.extend([(ct, v) for v in vals])
df = pd.DataFrame(rows, columns=["cell_type", "signature_expr"])

sns.set_theme(context="talk", style="white")
fig, ax = plt.subplots(figsize=(8, 5))

# Violin plot with quartile lines
sns.violinplot(
    data=df,
    x="cell_type",
    y="signature_expr",
    hue="cell_type",
    inner="quartile",     # ✅ quartile lines inside violins
    cut=0,
    density_norm="width", # ✅ modern replacement for scale="width"
    linewidth=1,
    palette="viridis",
    legend=False,
    ax=ax
)

# Overlay points (these look like “rings” outside quartiles)
sns.stripplot(
    data=df,
    x="cell_type",
    y="signature_expr",
    order=cell_types,
    color="black",
    size=3,
    jitter=0.15,
    alpha=0.4,
    ax=ax
)

ax.set_title("scRNA-seq huBCC", fontsize=16)
ax.set_ylabel("Expression level (normalized counts +1)")
ax.set_xlabel("")
sns.despine(offset=5, trim=True)
plt.tight_layout()
plt.show()