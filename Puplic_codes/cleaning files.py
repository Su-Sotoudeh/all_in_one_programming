'''import os
import pandas as pd

folder_path = 'C:\\Genomics\\results\\Genes_chardaff_parameters'
files = [f for f in os.listdir(folder_path) if f.endswith('.csv')]

df_list = []

for file in files:
    path = os.path.join(folder_path, file)

    # Load CSV
    df = pd.read_csv(path,low_memory=False)

    # Drop 'Unnamed' columns
    df.drop(columns=[col for col in df.columns if "Unnamed" in col], inplace=True)

    # Drop columns containing 'E_num'
    df.drop(columns=[col for col in df.columns if "E_num" in col], inplace=True)

    df_list.append(df)
    #print(f"{file}: {df.shape}")
    print(df.head(1).to_string())

# Concatenate all DataFrames
whole_csv = pd.concat(df_list, ignore_index=True)

# Save the concatenated CSV
whole_csv.to_csv(os.path.join(folder_path, "Whole_csv_file.csv"), index=False)
print(f"Whole CSV saved: {whole_csv.shape}")'''




# Save the concatenated CSV
import pandas as pd
df = pd.read_csv(r"C:\Genomics\\results\Genes_chardaff_parameters\Whole_csv_file.csv",low_memory=False)


# If Share_P% == 0 AND Share_C% == 0 → set Share_P2C = 0
df.loc[(df['CDS_P%'] == 0) & (df['CDS_C%'] == 0), 'CDS_P2C'] = 0
df.loc[(df['Share_P%'] == 0) & (df['Share_C%'] == 0), 'Share_P2C'] = 0
df.loc[(df['Int_P%'] == 0) & (df['Int_C%'] == 0), 'Int_P2C'] = 0
df.loc[(df['UTR_P%'] == 0) & (df['UTR_C%'] == 0), 'UTR_P2C'] = 0
maxraj_portion =   df['Int_ano'] +  df['UTR_ano']  + df['CDS_ano'] + df['share_ano']


Y = (df["Int_P_A2T"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_P_A2T"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_P_A2T"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_P_A2T"]*(df['share_ano']/maxraj_portion))
X = (df["Int_P2C"]   * (df['Int_ano']/maxraj_portion)) + (df["UTR_P2C"]  *(df['UTR_ano']/maxraj_portion)) +(df["CDS_P2C"]  *(df['CDS_ano']/maxraj_portion)) + (df["Share_P2C"]  *(df['share_ano']/maxraj_portion))




L = df["Length(M)"]
#  (df['CDS_P_A2T'] * (df['CDS_P%']/(df['CDS_P%']+df['CDS_C%']))) + (df['CDS_C_A2T'] * (df['CDS_C%']/(df['CDS_P%']+df['CDS_C%'])))
M_W_real = (df["Int_A2T_ano"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_A2T_ano"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_A2T_ano"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_A2T_ano"]*(df['share_ano']/maxraj_portion))

df['A2T_P'] = Y
df['P2C_W'] = X #annotated
df['A2T_annot'] = M_W_real
df.to_csv(r'C:\Genomics\\results\Genes_chardaff_parameters\Whole_csv_file.csv',index=False)


