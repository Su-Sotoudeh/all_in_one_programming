import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

eu_pro = pd.read_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\whole_datafram\Whole_datafram_2026.csv',low_memory=False)
eu_pro['Molecule_type'] = 'dsDNA'

Viruses = pd.read_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\Viruses\Viruses_Refrence_Datafram.csv',low_memory=False)

Plasmid = pd.read_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\Plasmides\Plasmides_refrence_Int_UTR_CDS_2026-02-05_23-22-23\Plasmides_2026-02-05_23-22-23_add_columns.csv',low_memory=False)


W_df = pd.concat([eu_pro,Viruses,Plasmid],axis=0,ignore_index=True)
#W_df = W_df.drop(columns='Assembly')
W_df.rename(columns={'Molecule_type':'Genome_Type'},inplace=True)


conditions = [
    (W_df['annot_P_A'] != 0) & (W_df['annot_C_A'] == 0),
    (W_df['annot_P_A'] == 0) & (W_df['annot_C_A'] != 0),
    (W_df['annot_P_A'] != 0) & (W_df['annot_C_A'] != 0)
]

choices = ['+', '-', '+/-']

W_df['Coding'] = np.select(conditions, choices, default='NA')


W_df = W_df[(W_df['annot_P_A'] != 0) | (W_df['annot_C_A'] != 0)]


print(W_df.shape)
print(W_df.head().to_string())


W_df.to_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\whole_datafram\Whole_datafram_2026_Eu_Pro_Virus_Plasmide.csv',index=False)