import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


sequence_file_path = r"/media/su/Bioinf/Genomes_2025_3/Genomes_2025_3_7/Chlroplast/ncbi_dataset/data/genome.fna"
df = pd.read_csv(r"/media/su/Bioinf/Genomes_2025_3/Genomes_2025_3_7/Chlroplast/ncbi_dataset/data/organelle_gene_report.csv",nrows=100)
print(df.head())
plus_list = []
minis_list = []
p_to_c = []
counter = 0
P2C_length_list = []
p2c_length = []
As = []
Ts = []
Cs = []
Gs = []
from Bio import SeqIO
for record in SeqIO.parse(sequence_file_path,"fasta"):
    genes_sequence_plus = ''
    genes_sequence_minis = ''
    id = record.id
    sequence = record.seq

    if counter == 100:
        break
    if id.startswith("NC"):
        counter+=1
        #print(record.id)
        print(counter)
        one_mt_df = df[df["accession"]==id]
        plus = 0
        minis = 0
        P_length = 0
        C_length = 0
        for index, row in one_mt_df.iterrows():
            #print(index, row["orientation"])
            start = row["start"]-1
            end = row["end"]
            gene_sequecne = sequence[start:end]

            length = int(end) - int(start)
            persent = length
            if row["orientation"] == "plus":
                genes_sequence_plus += gene_sequecne
                plus += 1
                P_length+= length
            else:
                genes_sequence_minis += gene_sequecne
                minis +=1
                C_length += length
        if plus == 0 and minis == 0 :
            continue
        #print('+',plus,'-',minis)
        P2C_length_list.append(P_length/(P_length+C_length)*100)
        plus_list.append(plus)
        minis_list.append(minis)
        p_to_c.append(plus/(plus+minis))
        A = genes_sequence_plus.count('A')/P_length*100
        T = genes_sequence_plus.count("T")/P_length*100
        C = genes_sequence_plus.count("C")/P_length*100
        G = genes_sequence_plus.count("G")/P_length*100
        As.append(A)
        Ts.append(T)
        Cs.append(C)
        Gs.append(G)
#x = np.arange(0,len(P2C_length_list))
#plt.scatter(x,plus_list)
#plt.scatter(x,minis_list)








import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Set global font to Times New Roman
#plt.rcParams["font.family"] = "Times New Roman"

# Example lists (replace with your data)
list1 = P2C_length_list
list2 = As
list3 = Ts
list4 = Cs
list5 = Gs

# Build dataframe
df = pd.DataFrame({
    "value": list1 + list2 + list3 + list4 + list5,
    "group": (["P/M"]*len(list1) +
              ["A"]*len(list2) +
              ["T"]*len(list3) +
              ["C"]*len(list4) +
              ["G"]*len(list5))
})

# Compute stats
stats = df.groupby("group")["value"].agg(["mean", "std", "var"])
print(stats)
print(df)
plt.figure(figsize=(8, 5))

# Custom palette
palette = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]

sns.violinplot(x="group", y="value", data=df, inner="box", palette=palette)

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
plt.ylim(0,120)
plt.tight_layout()


plt.savefig("mtgenome_violon.pdf", dpi=600, bbox_inches="tight")
plt.savefig("mtgenome_violon.png", dpi=600, bbox_inches="tight")

plt.show()












'''import pandas as pd

# Load JSONL
df = pd.read_json(r"/media/su/Bioinf/Genomes_2025_3/Genomes_2025_3_7/Chlroplast/ncbi_dataset/data/organelle_gene_report.jsonl", lines=True)

records = []
for row in df.to_dict(orient="records"):
    ann = row["annotations"][0]  # first annotation
    loc = ann["genomicLocations"][0]  # first location

    record = {

        "type": row["type"],
        "start": int(loc["genomicRange"]["begin"]),
        "end": int(loc["genomicRange"]["end"]),
        "orientation": loc["genomicRange"]["orientation"],
        "accession": loc["genomicAccessionVersion"],"taxname":row["taxname"],
        "taxId":row["taxId"],"assemblyAccession":ann["assemblyAccession"]
    }
    records.append(record)

flat_df = pd.DataFrame(records)
flat_df.to_csv(r"/media/su/Bioinf/Genomes_2025_3/Genomes_2025_3_7/Chlroplast/ncbi_dataset/data/organelle_gene_report.csv",index=False)

#print(df.head(30).to_string())'''

