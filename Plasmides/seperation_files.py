import os
from Bio import SeqIO
import pandas as pd

# Input files
fna_files = ["C:\Genomics\Genomes\Plasmid\plasmid.1.1.genomic.fna", "C:\Genomics\Genomes\Plasmid\plasmid.1.2.genomic.fna"]
gff_file = "C:\Genomics\Genomes\Plasmid\plasmid.1.genomic.gff"

# Output main folder
main_folder = "plasmides"
os.makedirs(main_folder, exist_ok=True)

# 1️⃣ Load GFF
gff_df = pd.read_csv(gff_file, sep="\t", comment="#", header=None,
                     names=["seqid","source","type","start","end","score","strand","phase","attributes"])

# 2️⃣ Load all sequences into a dict
seq_dict = {}
for fna in fna_files:
    for record in SeqIO.parse(fna, "fasta"):
        seq_dict[record.id] = record

# 3️⃣ Loop over all plasmid IDs in GFF
plasmid_ids = gff_df["seqid"].unique()

for pid in plasmid_ids:
    # Create folder for each plasmid
    folder_path = os.path.join(main_folder, pid)
    os.makedirs(folder_path, exist_ok=True)

    # Extract sequence
    if pid in seq_dict:
        SeqIO.write(seq_dict[pid], os.path.join(folder_path, f"{pid}.fna"), "fasta")
    else:
        print(f"Warning: sequence {pid} not found in FNA files!")

    # Extract GFF rows
    plasmid_gff = gff_df[gff_df["seqid"] == pid]
    plasmid_gff.to_csv(os.path.join(folder_path, f"{pid}.gff"),
                       sep="\t", header=False, index=False)
