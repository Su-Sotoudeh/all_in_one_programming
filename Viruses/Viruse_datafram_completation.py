
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


'''path = r'C:\Genomics\Chargaff\coding\all_in_one_programming\Viruses\Viruses_refrence_Int_UTR_CDS_2026-01-26_16-51-15\Viruses_2026-01-26_16-51-15.csv'
df_v =  pd.read_csv(path,low_memory=False)
df_v['Group'] = 'Viruses'
df = df_v[df_v['Wh_ano']>0]

#df = df_v.drop(columns=[	'Int_P_E_num','Int_C_E_num','UTR_P_E_num','UTR_C_E_num','CDS_P_E_num','CDS_C_E_num','Share_P_E_num','Share_C_E_num','Unanot_P_E_num','Unanot_C_E_num' ])


# Find rows where both conditions are True
mask = (df['UTR_P%'] == 0) & (df['UTR_C%'] == 0)
mask2 = (df['Share_P%'] == 0) & (df['Share_C%'] == 0)
mask3 = (df['Int_P%'] == 0) & (df['Int_C%'] == 0)
mask4 = (df['CDS_P%'] == 0) & (df['CDS_C%'] == 0)
# Set 'UTR_P2C' to 0 only for those rows
df.loc[mask, 'UTR_P2C'] = 0
df.loc[mask2, 'Share_P2C'] = 0
df.loc[mask3, 'Int_P2C'] = 0
df.loc[mask4, 'CDS_P2C'] = 0


maxraj_portion =   df['Int_ano'] +  df['UTR_ano']  + df['CDS_ano'] + df['share_ano']

A2T_P = (df["Int_P_A2T"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_P_A2T"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_P_A2T"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_P_A2T"]*(df['share_ano']/maxraj_portion))
A2T_C = (df["Int_C_A2T"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_C_A2T"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_C_A2T"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_C_A2T"]*(df['share_ano']/maxraj_portion))

G2C_P = (df["Int_P_G2C"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_P_G2C"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_P_G2C"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_P_G2C"]*(df['share_ano']/maxraj_portion))
G2C_C = (df["Int_C_G2C"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_C_G2C"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_C_G2C"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_C_G2C"]*(df['share_ano']/maxraj_portion))


P2C_W = (df["Int_P2C"]   * (df['Int_ano']/maxraj_portion)) + (df["UTR_P2C"]  *(df['UTR_ano']/maxraj_portion)) +(df["CDS_P2C"]  *(df['CDS_ano']/maxraj_portion)) + (df["Share_P2C"]  *(df['share_ano']/maxraj_portion))

#  (df['CDS_P_A2T'] * (df['CDS_P%']/(df['CDS_P%']+df['CDS_C%']))) + (df['CDS_C_A2T'] * (df['CDS_C%']/(df['CDS_P%']+df['CDS_C%'])))

A2T_annotation = (df["Int_A2T_ano"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_A2T_ano"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_A2T_ano"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_A2T_ano"]*(df['share_ano']/maxraj_portion))
G2C_annotation = (df["Int_G2C_ano"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_G2C_ano"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_G2C_ano"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_G2C_ano"]*(df['share_ano']/maxraj_portion))

Ap_annotation = (df["Int_P_A_ano"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_P_A_ano"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_P_A_ano"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_P_A_ano"]*(df['share_ano']/maxraj_portion))
Tp_annotation = (df["Int_P_T_ano"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_P_T_ano"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_P_T_ano"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_P_T_ano"]*(df['share_ano']/maxraj_portion))
Cp_annotation = (df["Int_P_C_ano"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_P_C_ano"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_P_C_ano"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_P_C_ano"]*(df['share_ano']/maxraj_portion))
Gp_annotation = (df["Int_P_G_ano"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_P_G_ano"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_P_G_ano"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_P_G_ano"]*(df['share_ano']/maxraj_portion))

Ac_annotation = (df["Int_C_A_ano"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_C_A_ano"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_C_A_ano"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_C_A_ano"]*(df['share_ano']/maxraj_portion))
Tc_annotation = (df["Int_C_T_ano"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_C_T_ano"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_C_T_ano"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_C_T_ano"]*(df['share_ano']/maxraj_portion))
Cc_annotation = (df["Int_C_C_ano"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_C_C_ano"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_C_C_ano"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_C_C_ano"]*(df['share_ano']/maxraj_portion))
Gc_annotation = (df["Int_C_G_ano"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_C_G_ano"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_C_G_ano"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_C_G_ano"]*(df['share_ano']/maxraj_portion))


df["annot_P_A"] = Ap_annotation
df["annot_P_T"] = Tp_annotation
df["annot_P_C"] = Cp_annotation
df["annot_P_G"] = Gp_annotation

df["annot_C_A"] = Ac_annotation
df["annot_C_T"] = Tc_annotation
df["annot_C_C"] = Cc_annotation
df["annot_C_G"] = Gc_annotation

df['A2T_P'] = A2T_P
df['A2T_C'] = A2T_C

df['G2C_P'] = G2C_P
df['G2C_C'] = G2C_C


df['A2T_annot'] = A2T_annotation
df['G2C_annot'] = G2C_annotation

df['P2C_W'] = P2C_W #annotated


df['ATp_Skew'] = (Ap_annotation - Tc_annotation)/(Ap_annotation + Tp_annotation)
df['GCp_Skew'] = (Gp_annotation - Cp_annotation)/(Gp_annotation + Cp_annotation)
df['ATc_Skew'] = (Ac_annotation - Tc_annotation)/(Ap_annotation + Tc_annotation)
df['GCc_Skew'] = (Gc_annotation - Cc_annotation)/(Gc_annotation + Cc_annotation)
df['ATpc_Skew'] = ((Ap_annotation + Ac_annotation) - (Tp_annotation + Tc_annotation))/((Ap_annotation + Ac_annotation) + (Tp_annotation + Tc_annotation))
df['GCpc_Skew'] = ((Gp_annotation + Gc_annotation) - (Cp_annotation + Cc_annotation))/((Gp_annotation + Gc_annotation) + (Cp_annotation + Cc_annotation))

df['A2T_P-A2T_C'] = df['A2T_P'] - df['A2T_C']
df['G2C_P-G2C_C'] = df['G2C_P'] - df['G2C_C']

print(df.isnull().sum())

df.to_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\Viruses\Viruses_dstafram_2026-01-26_added_columns.csv'
          ,index=False)
df = pd.read_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\Viruses\Viruses_dstafram_2026-01-26_added_columns.csv')

#df =df[(df['Chr_ID'].str.startswith('NC'))&(df['Wh_ano']>10)]
df = df.sort_values('A2T_P')

A2T_P = df['A2T_P']
x = np.arange(len(A2T_P))

plt.scatter(x,A2T_P)
plt.show()'''


#----------------------------------------------------------------
#                       Masking for Refrence genomes
#-------------------------------------------------------------------



import pandas as pd

# Load main dataframe
path = r'C:\Genomics\Chargaff\coding\all_in_one_programming\Viruses\Viruses_refrence_Int_UTR_CDS_2026-01-26_16-51-15\Viruses_2026-01-26_16-51-15.csv'
df_v = pd.read_csv(path, low_memory=False)

# Drop unwanted column
df_v = df_v.drop(columns=['Genome'], errors='ignore')

# Load reference dataframe
ref_df = pd.read_csv(r"C:\Genomics\Genomes\unzip_viruse_all_complete\unzip_viruse_all_complete\Viruse_file_information\Refrence_viruses_list.csv")

# Keep only needed columns
ref_df = ref_df[["Assembly", "Molecule_type"]]

# Merge main df with reference df
merged_df = df_v.merge(
    ref_df,
    left_on="Name",
    right_on="Assembly",
    how="inner"     # keeps only rows that exist in reference
)

# Drop the duplicate key column if you don't want it
merged_df = merged_df.drop(columns=["Assembly"])

# Save or inspect
print(merged_df.head().to_string())



print(merged_df.shape)
df = merged_df

merged_df.to_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\Viruses\Viruses_Refrence_Datafram.csv',index=False)

