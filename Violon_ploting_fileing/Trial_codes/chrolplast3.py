#this is my code that need to be sutabile
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns




sequence_file_path = r"/media/su/Bioinf/Genomes_2025_3/Genomes_2025_3_7/Chlroplast/ncbi_dataset/data/genome.fna"
df = pd.read_csv(r"/media/su/Bioinf/Genomes_2025_3/Genomes_2025_3_7/Chlroplast/ncbi_dataset/data/organelle_gene_report.csv",nrows=10000)
print(df.head())


# Define Okabe-Ito color palette (colorblind-friendly)
okabe_ito = ['#9467bd', '#2ca02c','#1f77b4', '#ff7f0e','#d62728', ]  # Blue, Orange, Green, Red, Purple
result_p = {'P_poerion': None, 'Ap': None, 'Tp': None, 'Cp': None, 'Gp': None}
result_m = {'M_portion': None, 'Ac': None, 'Tc': None, 'Cc': None, 'Gc': None}
P_portion_s = []
Aps = []
Tps = []
Cps = []
Gps = []
M_portion_s = []
Ams = []
Tms = []
Cms = []
Gms = []

# Function to process a single record (for parallelization)
counter = 0
from Bio import SeqIO
for record in SeqIO.parse(sequence_file_path,"fasta"):
    counter += 1
    print(counter)
    #if counter == 100000:
    #    break
    genes_sequence_plus = ''
    genes_sequence_minis = ''
    id = record.id
    sequence = record.seq


    if id.startswith("NC"):
        one_mt_df = df[df["accession"] == id]
        plus = 0
        minis = 0
        P_length = 0
        M_length = 0
        for _, row in one_mt_df.iterrows():
            start = row["start"] - 1
            end = row["end"]
            gene_sequence = sequence[start:end]
            length = end - start
            if row["orientation"] == "plus":
                genes_sequence_plus += gene_sequence
                plus += 1
                P_length += length
            else:
                genes_sequence_minis += gene_sequence
                minis += 1
                M_length += length
        if plus == 0 and minis == 0:
            continue
        P_portion = len(genes_sequence_plus) / len(genes_sequence_plus + genes_sequence_minis) * 100 if(len(genes_sequence_plus) + len(genes_sequence_minis)) > 0 else 0
        Ap = genes_sequence_plus.count('A') / P_length * 100 if P_length > 0 else 0
        Tp = genes_sequence_plus.count('T') / P_length * 100 if P_length > 0 else 0
        Cp = genes_sequence_plus.count('C') / P_length * 100 if P_length > 0 else 0
        Gp = genes_sequence_plus.count('G') / P_length * 100 if P_length > 0 else 0

        M_portion = len(genes_sequence_minis) / len(genes_sequence_plus + genes_sequence_minis) * 100 if(len(genes_sequence_plus) + len(genes_sequence_minis)) > 0 else 0
        Am = genes_sequence_minis.count('A') / P_length * 100 if P_length > 0 else 0
        Tm = genes_sequence_minis.count('T') / P_length * 100 if P_length > 0 else 0
        Cm = genes_sequence_minis.count('C') / P_length * 100 if P_length > 0 else 0
        Gm = genes_sequence_minis.count('G') / P_length * 100 if P_length > 0 else 0
        #results = {'P_portion': P_portion, 'Ap': Ap, 'Tp': Tp, 'Cp': Cp, 'Gp': Gp,'M_portion': M_portion, 'Am': Am, 'Tm': Tm, 'Cm': Cm, 'Gm': Gm}
        P_portion_s.append(P_portion)
        Aps.append(Ap)
        Tps.append(Tp)
        Cps.append(Cp)
        Gps.append(Gp)
        M_portion_s.append(M_portion)
        Ams.append(Am)
        Tms.append(Tm)
        Cms.append(Cm)
        Gms.append(Gm)



list1 = P_portion_s
list2 = Aps
list3 = Tps
list4 = Cps
list5 = Gps
list6 = M_portion_s
list7 = Ams
list8 = Tms
list9 = Cms
list10 =Gms






# Build dataframe
df = pd.DataFrame({
    "value": list1 + list2 + list3 + list4 + list5+list6+list7+list8+list9+list10 ,
    "group": (["P"]*len(list1) +
              ["Ap"]*len(list2)   +
              ["Tp"]*len(list3)   +
              ["Cp"]*len(list4)   +
              ["Gp"]*len(list5)   +
              ["M"]*len(list1) +
              ["Am"]*len(list2)   +
              ["Tm"]*len(list3)   +
              ["Cm"]*len(list4)   +
              ["Gm"]*len(list5)

                                    )
})

# Compute stats
stats = df.groupby("group")["value"].agg(["mean", "std", "var"])
print(stats)

plt.figure(figsize=(10, 5))

# Custom palette
palette = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd","#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]

sns.violinplot(x="group", y="value", data=df, inner="box", palette=okabe_ito)

# Add mean and std as text above each violin
for i, group in enumerate(df["group"].unique()):
    m = stats.loc[group, "mean"]
    s = stats.loc[group, "std"]
    ymax = df[df["group"] == group]["value"].max()

    plt.text(i, ymax*1.4, f"Mean={m:.2f}\nSD={s:.2f}",
             horizontalalignment='center', color="black", weight="bold")

plt.title("Violin Plot with Mean and SD", fontsize=14)
plt.xlabel("Groups", fontsize=12)
plt.ylabel("Values", fontsize=12)
plt.ylim(0,100)
plt.tight_layout()






plt.savefig("mtgenome_violon.pdf", dpi=600, bbox_inches="tight")
plt.savefig("mtgenome_violon.png", dpi=600, bbox_inches="tight")

plt.show()