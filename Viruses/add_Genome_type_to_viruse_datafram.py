import pandas as pd

# 1. Load both dataframes
main_path = r'C:\Genomics\Chargaff\coding\all_in_one_programming\Viruses\Viruses_refrence_Int_UTR_CDS_2026-01-26_16-51-15\Viruses_2026-01-26_16-51-15.csv'
info_path = r'C:\Genomics\Chargaff\coding\all_in_one_programming\Viruses\Viruse_file_information\Viruses_Genome_type.csv'
#df_info = pd.read_csv(info_path, low_memory=False)
#print(df_info.isnull().sum())


df_main = pd.read_csv(main_path, low_memory=False)
print(df_main.isnull().sum())
if 'Assembly Accession' in df_main.columns:
    df_main = df_main.drop(columns=['Assembly Accession','Assembly Accession_x','Assembly Accession_y'])
if 'Genome' in df_main.columns:
    df_main = df_main.drop(columns=['Genome','Genome_x','Genome_y'])


df_info = pd.read_csv(info_path, low_memory=False)
print(df_main.head().to_string())
# 2. Select ONLY the 'Assembly Accession' and 'Genome' columns
df_info_subset = df_info[['Assembly Accession', 'Genome']]

# 3. Merge them
# Matching 'Name' from main with 'Assembly Accession' from info
df_merged = pd.merge(
    df_main,
    df_info_subset,
    left_on='Name',
    right_on='Assembly Accession',
    how='left'
)

df_merged = df_merged.drop(columns=['Assembly Accession'])

# 5. Check the result
print(f"Original shape: {df_main.shape}")
print(f"Merged shape: {df_merged.shape}")
print("New columns:", df_merged.columns)






print(df_merged.head())
print(df_info.isnull().sum())
df_merged.to_csv(main_path, index=False)

print("File saved successfully.")