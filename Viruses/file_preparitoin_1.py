# here we have fixing Viruses_information_NCBI_whole.csv file
#we shoud add it phylogen innformmation
# we used taxon ID to optain phyo nemes infrommation
# final as output we have ncbi_with_taxonomy_parallel.csv file





import pandas as pd
from Bio import Entrez
from multiprocessing import Pool, cpu_count
import time

# REQUIRED by NCBI
Entrez.email = "reza.sotoudeh.van@gmail.com"

# Load your NCBI CSV
df = pd.read_csv(r"D:\PhD\genomes\unzip_viruse_all_complete\Viruses_information_NCBI_whole.csv")
print(df.shape)
def fetch_taxonomy(taxid):
    """Fetch taxonomy lineage for a given TaxID from NCBI Taxonomy database."""
    try:
        handle = Entrez.efetch(db="taxonomy", id=str(taxid), retmode="xml")
        records = Entrez.read(handle)
        lineage_info = records[0]['LineageEx']
        lineage_dict = {rank['Rank']: rank['ScientificName'] for rank in lineage_info}
        lineage_dict['Species'] = records[0]['ScientificName']
        lineage_dict['TaxID'] = taxid
        print(lineage_dict)
        return lineage_dict
    except Exception as e:
        print(f"❌ Error fetching TaxID {taxid}: {e}")
        return {"TaxID": taxid}

def parallel_fetch(taxids, workers=4):
    """Fetch taxonomy in parallel with multiple processes."""
    with Pool(processes=workers) as pool:
        results = pool.map(fetch_taxonomy, taxids)
    return results

if __name__ == "__main__":
    taxid_list = df["Organism Taxonomic ID"].dropna().unique().tolist()

    print(f"🔎 Fetching taxonomy for {len(taxid_list)} TaxIDs...")

    # Use min(cpu_count(), 4) to avoid hammering NCBI
    results = parallel_fetch(taxid_list, workers=min(cpu_count(), 4))

    # Convert list of dicts → DataFrame
    taxonomy_df = pd.DataFrame(results)

    # Merge back into original DataFrame
    merged = pd.merge(df, taxonomy_df, left_on="Organism Taxonomic ID", right_on="TaxID", how="left")

    # Save output
    merged.to_csv("ncbi_with_taxonomy_parallel.csv", index=False)
    print("✅ Done! Saved to ncbi_with_taxonomy_parallel.csv")