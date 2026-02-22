import numpy as np
import pandas as pd

# Load your data
my_datafram = pd.read_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\whole_datafram\Whole_datafram_2026_Eu_Pro_Virus_Plasmide.csv', low_memory=False)
print(my_datafram.head().to_string())
# List to store data for the table
data_for_table = []

groups = my_datafram['Group'].unique()
print(groups)
for group in groups:

    # --- SPECIAL CASE FOR VIRUSES ---
    if group == "Viruses":
        coding_types = ["+", "-", "+/-"]
        subgroups = [(group + "_" + c, my_datafram[(my_datafram['Group'] == group) & (my_datafram['Coding'] == c)])
                     for c in coding_types]
    else:
        subgroups = [(group, my_datafram[my_datafram['Group'] == group])]

    # --- PROCESS EACH SUBGROUP ---
    for subgroup_name, df in subgroups:

        if df.empty:
            continue

        # Calculate Length
        L = sum(df['Length(M)'])

        # Percentages
        int_val = (sum(df['Int_ano']/100 * df['Length(M)']) / L) * 100
        utr_val = (sum(df['UTR_ano']/100 * df['Length(M)']) / L) * 100
        cds_val = (sum(df['CDS_ano']/100 * df['Length(M)']) / L) * 100
        share_val = (sum(df['share_ano']/100 * df['Length(M)']) / L) * 100
        Un_ano_val = (sum(df['Un_ano']/100 * df['Length(M)']) / L) * 100

        # Portions
        P_annotaion_portion = sum(
            (df['Int_P%']* df['Length(M)']) +
            (df['UTR_P%']*df['Length(M)']) +
            (df['CDS_P%']*df['Length(M)']) +
            (df['Share_P%']*df['Length(M)'])
        ) / 100

        C_annotaion_portion = sum(
            (df['Int_C%']* df['Length(M)']) +
            (df['UTR_C%']* df['Length(M)']) +
            (df['CDS_C%']* df['Length(M)']) +
            (df['Share_C%']* df['Length(M)'])
        ) / 100

        Un_ano_portion = sum(df['Un_ano']/100 * df['Length(M)'])

        organism_number = len(df['Name'].unique())
        chr_number = len(df['Chr_ID'].unique())

        # Sum/L ratio
        sum_L_ratio = (Un_ano_portion + P_annotaion_portion + C_annotaion_portion) / L

        # Store results
        row_data = {
            'Group': subgroup_name,
            'Organ_Num': organism_number,
            'Chr_Num': chr_number,
            'Length (Bbp)': round(L/1000, 3),
            '(+)Annot (Bbp)': round(P_annotaion_portion/1000, 3),
            '(-)Annot (Bbp)': round(C_annotaion_portion/1000, 3),
            'Unannot (Bbp)': round(Un_ano_portion/1000, 3),
            '(+/-)': round(P_annotaion_portion/C_annotaion_portion, 3) if C_annotaion_portion != 0 else None,
            'Int (%)': round(int_val, 3),
            'UTR (%)': round(utr_val, 3),
            'CDS (%)': round(cds_val, 3),
            'Share (%)': round(share_val, 3),
            'Unannot (%)': round(Un_ano_val, 3)
        }

        data_for_table.append(row_data)

        # Print text info
        print(f'--- Group: {subgroup_name} ---')
        print('organ number', organism_number)
        print('chr number', chr_number)
        print('whole sequence length is', round(L/1000, 3), 'Bbp')
        print('Unannot_portion', round(Un_ano_portion/1000, 3), 'Bbp')
        print('(+)annotation', round(P_annotaion_portion/1000, 3), 'Bbp')
        print('(-)annotation', round(C_annotaion_portion/1000, 3), 'Bbp')
        print('(+/-)', round(P_annotaion_portion/C_annotaion_portion, 3) if C_annotaion_portion != 0 else "NA")
        print('Int', round(int_val, 3), '%')
        print('UTR', round(utr_val, 3), '%')
        print('CDS', round(cds_val, 3), '%')
        print('Share', round(share_val, 3), '%')
        print('Un annot', round(Un_ano_val, 3), '%')
        print('-----------------------')

# --- CREATE THE TABLE ---
results_df = pd.DataFrame(data_for_table)

print("\nSummary Table:")
print(results_df.to_string())

results_df.to_csv('Group_Summary_Table_v2.csv', index=False)

print("\nTable saved to 'Group_Summary_Table_v2.csv'")


import numpy as np
import pandas as pd

# Load your data
my_datafram = pd.read_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\whole_datafram\Whole_datafram_2026_Eu_Pro_Virus_Plasmide.csv', low_memory=False)

# List to store data for the table
data_for_table = []

groups = my_datafram['Group'].unique()

for group in groups:

    # --- SPECIAL CASE FOR VIRUSES ---
    if group == "Viruses":

        genome_types = ["DNA", "RNA"]
        coding_types = ["+", "-", "+/-"]

        subgroups = []

        for gtype in genome_types:
            for ctype in coding_types:
                df_sub = my_datafram[
                    (my_datafram['Group'] == group) &
                    (my_datafram['Genome_Type'].str.contains(gtype, case=False, na=False)) &
                    (my_datafram['Coding'] == ctype)
                    ]

                subgroup_name = f"{group}_{gtype}_{ctype}"
                subgroups.append((subgroup_name, df_sub))

    else:
        # Normal groups
        df_sub = my_datafram[my_datafram['Group'] == group]
        subgroups = [(group, df_sub)]

    # --- PROCESS EACH SUBGROUP ---
    for subgroup_name, df in subgroups:

        if df.empty:
            continue

        # Calculate Length
        L = sum(df['Length(M)'])

        # Percentages
        int_val = (sum(df['Int_ano']/100 * df['Length(M)']) / L) * 100
        utr_val = (sum(df['UTR_ano']/100 * df['Length(M)']) / L) * 100
        cds_val = (sum(df['CDS_ano']/100 * df['Length(M)']) / L) * 100
        share_val = (sum(df['share_ano']/100 * df['Length(M)']) / L) * 100
        Un_ano_val = (sum(df['Un_ano']/100 * df['Length(M)']) / L) * 100

        # Portions
        P_annotaion_portion = sum(
            (df['Int_P%']* df['Length(M)']) +
            (df['UTR_P%']*df['Length(M)']) +
            (df['CDS_P%']*df['Length(M)']) +
            (df['Share_P%']*df['Length(M)'])
        ) / 100

        C_annotaion_portion = sum(
            (df['Int_C%']* df['Length(M)']) +
            (df['UTR_C%']* df['Length(M)']) +
            (df['CDS_C%']* df['Length(M)']) +
            (df['Share_C%']* df['Length(M)'])
        ) / 100

        Un_ano_portion = sum(df['Un_ano']/100 * df['Length(M)'])

        organism_number = len(df['Name'].unique())
        chr_number = len(df['Chr_ID'].unique())

        # Store results
        d = 6
        row_data = {
            'Group': subgroup_name,
            'Organ_Num': organism_number,
            'Chr_Num': chr_number,
            'Length (Bbp)': round(L/1000, d),
            '(+)Annot (Bbp)': round(P_annotaion_portion/1000, d),
            '(-)Annot (Bbp)': round(C_annotaion_portion/1000, d),
            'Unannot (Bbp)': round(Un_ano_portion/1000, d),
            '(+/-)': round(P_annotaion_portion/C_annotaion_portion, d) if C_annotaion_portion != 0 else None,
            'Int (%)': round(int_val, d),
            'UTR (%)': round(utr_val, d),
            'CDS (%)': round(cds_val, d),
            'Share (%)': round(share_val, d),
            'Unannot (%)': round(Un_ano_val, d)
        }

        data_for_table.append(row_data)

        # Print text info
        print(f'--- Group: {subgroup_name} ---')
        print('organ number', organism_number)
        print('chr number', chr_number)
        print('whole sequence length is', round(L/1000, d), 'Bbp')
        print('Unannot_portion', round(Un_ano_portion/1000, d), 'Bbp')
        print('(+)annotation', round(P_annotaion_portion/1000, d), 'Bbp')
        print('(-)annotation', round(C_annotaion_portion/1000, d), 'Bbp')
        print('(+/-)', row_data['(+/-)'])
        print('Int', round(int_val, d), '%')
        print('UTR', round(utr_val, d), '%')
        print('CDS', round(cds_val, d), '%')
        print('Share', round(share_val, d), '%')
        print('Un annot', round(Un_ano_val, d), '%')
        print('-----------------------')

# --- CREATE THE TABLE ---
results_df = pd.DataFrame(data_for_table)

print("\nSummary Table:")
print(results_df.to_string())

results_df.to_csv('Group_Summary_Table_v2.csv', index=False)

print("\nTable saved to 'Group_Summary_Table_v2.csv'")
