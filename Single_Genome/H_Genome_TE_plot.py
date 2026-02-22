import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Load your data
path = r'D:\PhD\genomes\Mus_19\NC_000085.7.csv'
#path = '/media/su/Bioinf/output/Dolphin/NC_047034.1.csv'
df = pd.read_csv(path)
print(df)
# Extract data
labels = ['CDS', 'UTR', 'Intron', 'LINE', 'SINE', 'LTR', 'DNA']
seri_plus = df['whole'][18:25].values
seri_minis = df['whole'][25:32].values
# Normalize to percentage
seri_plus_pct = seri_plus / seri_plus.sum() * 100
seri_minis_pct = seri_minis / seri_minis.sum() * 100

# Set up bar plot
x = np.arange(len(labels))
width = 0.35

# Set theme
sns.set_theme(style="whitegrid", font_scale=1.3)

# Choose a color palette
#colors = sns.color_palette("muted")  # Try "muted", "rocket", "colorblind", etc.
colors = sns.color_palette("colorblind")
# Plotting
fig, ax = plt.subplots(figsize=(10, 6))
bars1 = ax.bar(x - width/2, seri_plus_pct, width, label='+ annotated', color=colors[0])
bars2 = ax.bar(x + width/2, seri_minis_pct, width, label='- annotated', color=colors[1])

# Add value labels
for bar in bars1 + bars2:
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2.0, height + 0.5, f'{height:.1f}%', ha='center', va='bottom', fontsize=10)

# Labels and formatting
ax.set_ylabel('Percentage (%)')
ax.set_xlabel('Elements')
ax.set_title('Genome Element Percentages on + and - Strands')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()

# Clean up and show
sns.despine()
plt.tight_layout()
plt.savefig("genome_elements_percent_plot.png", dpi=600)
plt.show()
