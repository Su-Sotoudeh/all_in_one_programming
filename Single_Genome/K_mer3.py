# kmer_symmetry_pipeline.py
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

# ---------------- Load data ----------------
# First column = organism names (no header)
df = pd.read_csv("../112_60.csv", sep=",", header=0)   # adjust sep="," if CSV is comma-separated
df.rename(columns={df.columns[0]: "organism"}, inplace=True)

# Extract k-mer base names (before the _+M, _+RC, etc.)
kmers = sorted(set([c.split("_")[0] for c in df.columns if c != "organism"]))

# ---------------- Symmetry correlations ----------------
results = []
for k in kmers:
    plusM = df[f"{k}_+M"]
    minusRC = df[f"{k}_-RC"]
    minusM = df[f"{k}_-M"]
    plusRC = df[f"{k}_+RC"]

    # Correlations
    r1, _ = pearsonr(plusM, minusRC)   # (+M vs -RC)
    r2, _ = pearsonr(minusM, plusRC)   # (-M vs +RC)

    # Symmetry index
    S = ((abs(plusM - minusRC) + abs(minusM - plusRC)) /
         (plusM + minusM + plusRC + minusRC).replace(0, 1)).mean()

    results.append([k, r1, r2, S])

res_df = pd.DataFrame(results, columns=["kmer", "corr(+M,-RC)", "corr(-M,+RC)", "symmetry_index"])
print("\n=== Symmetry correlation summary ===")
print(res_df.sort_values("symmetry_index"))

# ---------------- Scatter plot (example) ----------------
k = "CA"
plt.figure(figsize=(6,5))
sns.scatterplot(x=df[f"{k}_+M"], y=df[f"{k}_-RC"], alpha=0.6)
sns.regplot(x=df[f"{k}_+M"], y=df[f"{k}_-RC"], scatter=False, color="red")
plt.xlabel(f"{k} +M")
plt.ylabel(f"{k} -RC")
plt.title(f"Scatter {k}: +M vs -RC")
plt.show()

# ---------------- Heatmap of Δ values ----------------
heat_data = []
for k in kmers:
    delta1 = df[f"{k}_+M"] - df[f"{k}_-RC"]
    delta2 = df[f"{k}_-M"] - df[f"{k}_+RC"]
    heat_data.append((k, delta1.mean(), delta2.mean()))

heat_df = pd.DataFrame(heat_data, columns=["kmer", "Δ(+M - -RC)", "Δ(-M - +RC)"])
plt.figure(figsize=(10,4))
sns.heatmap(heat_df.set_index("kmer").T, cmap="coolwarm", center=0, annot=True, fmt=".2f")
plt.title("Mean Δ deviations across organisms")
plt.show()

# ---------------- Symmetry index distribution ----------------
plt.figure(figsize=(8,4))
sns.barplot(x="kmer", y="symmetry_index", data=res_df, palette="viridis")
plt.xticks(rotation=45)
plt.ylabel("Symmetry Index (S)")
plt.title("Symmetry Index per k-mer (lower = more symmetric)")
plt.tight_layout()
plt.show()
