# here we have finall fix and compeletion on NCBI_ICTV_match_genus_genome_type.csv file
# some row lost its infromation or was empety they fixed with new fetching
# here we have viruses_with_taxonomy_fixed.csv file as output file





import pandas as pd
from Bio import Entrez

Entrez.email = "reza.sotoudeh.van@gmail.com"

# ---------- Step 1: Load your existing results ----------
df = pd.read_csv("ncbi_with_taxonomy_parallel.csv")

# Failed rows: where 'species' (or any rank column) is missing
failed = df[df["species"].isna()]["TaxID"].dropna().astype(str).unique().tolist()
print(f"❌ Found {len(failed)} TaxIDs missing taxonomy")

# ---------- Step 2: Function to fetch again ----------
def fetch_taxonomy_batch(taxid_list):
    try:
        handle = Entrez.efetch(db="taxonomy", id=",".join(taxid_list), retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        results = {}
        for rec in records:
            lineage = {item["Rank"]: item["ScientificName"] for item in rec["LineageEx"]}
            lineage["species"] = rec.get("ScientificName", "")
            lineage["TaxID"] = rec["TaxId"]
            results[rec["TaxId"]] = lineage
        return results
    except Exception as e:
        print("⚠️ Error batch:", e)
        return {tid: None for tid in taxid_list}

# ---------- Step 3: Retry in batches ----------
batch_size = 1000
retrieved = {}

for i in range(0, len(failed), batch_size):
    chunk = failed[i:i+batch_size]
    print(f"Fetching batch {i//batch_size+1} / {len(failed)//batch_size+1} ...")
    results = fetch_taxonomy_batch(chunk)
    retrieved.update(results)

# ---------- Step 4: Make dataframe of retries ----------
retry_df = pd.DataFrame.from_dict(retrieved, orient="index")
retry_df.reset_index(drop=True, inplace=True)

# ---------- Step 5: Merge into your old file ----------
for col in retry_df.columns:
    if col not in df.columns:
        df[col] = None

for _, row in retry_df.iterrows():
    taxid = str(row["TaxID"])
    for col in retry_df.columns:
        if pd.notna(row[col]):
            df.loc[df["TaxID"] == int(taxid), col] = row[col]

# Save updated file
df.to_csv("viruses_with_taxonomy_fixed.csv", index=False)
print("✅ Saved merged file as viruses_with_taxonomy_fixed.csv")