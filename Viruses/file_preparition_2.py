#here we use two file
#ncbi_with_taxonomy_parallel.csv from NCBI(and added phylum name by we
#viruses_genome_information.csv from ICTV
#ICTV have genopme Type and phlum information
#NCBI have Genome ID and phylum information
# we want to connect Genome ID to Genome type $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#here we by use they shered part(phylum) we connect Genomme ID and Genomye Type
# and as output we have NCBI_ICTV_match_genus_genome_type.csv


import pandas as pd

# Load files
ncbi_df = pd.read_csv(r"D:\Projects\all_in_one_programming\ncbi_with_taxonomy_parallel.csv")
ictv_df = pd.read_csv(r"D:\PhD\genomes\unzip_viruse_all_complete\viruses_genome_information.csv")


import pandas as pd

# ---------------- Load data ----------------
ncbi_df = pd.read_csv(
    r"D:\\Projects\\all_in_one_programming\\viruses_with_taxonomy_fixed.csv",
    low_memory=False,
    dtype=str
)
ictv_df = pd.read_csv(
    r"D:\PhD\genomes\\unzip_viruse_all_complete\\viruses_genome_information.csv",
    low_memory=False,
    dtype=str
)

# Fill missing values with empty strings
ncbi_df = ncbi_df.fillna("")
ictv_df = ictv_df.fillna("")

# ---------------- Prepare ICTV info ----------------
# We only need: genus → genome type + virus name
ictv_genus = ictv_df.groupby("Genus").first().reset_index()
ictv_genus = ictv_genus[["Genus", "Virus name(s)", "Genome"]]

# ---------------- Merge by genus ----------------
merged_df = pd.merge(
    ncbi_df,
    ictv_genus,
    left_on="genus",   # from NCBI file
    right_on="Genus",  # from ICTV file
    how="inner"        # keep only matched
)

# ---------------- Select final columns ----------------
final_df = merged_df[
    [
        "Assembly Accession",
        "Organism Taxonomic ID",
        "Organism Name",
        "realm",
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "Virus name(s)",
        "Genome"
    ]
]

# ---------------- Save result ----------------
output_path = r"D:\Projects\\all_in_one_programming\\NCBI_ICTV_match_genus_genome_type.csv"
final_df.to_csv(output_path, index=False)

print(f"✅ Finished! Matched {len(final_df)} viruses with genome type info.")
print(f"Saved to: {output_path}")