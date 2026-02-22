import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
my_datafram = pd.read_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\Viruses\Viruses_refrence_Int_UTR_CDS_2026-01-26_16-51-15\Viruses_2026-01-26_16-51-15.csv',low_memory=False)
#we fillter it  as refrence genomse

refrence = pd.read_csv(r"C:\Genomics\Genomes\unzip_viruse_all_complete\unzip_viruse_all_complete\Viruse_file_information\Refrence_viruses_list.csv")
print(refrence['Molecule_type'].unique())


main_df = my_datafram
ref_df  = refrence

# Keep only needed columns from reference
ref_df = ref_df[["Assembly", "Molecule_type"]]

# Merge: this automatically filters to RefSeq genomes
filtered = main_df.merge(
    ref_df,
    left_on="Name",
    right_on="Assembly",
    how="inner"
)
filtered["Genome"] = filtered["Molecule_type"]

my_datafram = filtered

print(my_datafram.head().to_string())
print(my_datafram.shape)
df = my_datafram

df = df[(df['annot_P_A'] != 0) & (df['annot_C_A'] == 0)]# condition = '+'
df = df[(df['annot_P_A'] == 0) & (df['annot_C_A'] != 0)]# condition = '-'
df = df[(df['annot_P_A'] != 0) & (df['annot_C_A'] != 0)]# condition = '+/-'

group = 'RNA'
group = 'DNA'
condition = '-'#

df = df[df["Genome"].str.contains(group, na=False)]





print(f"--- Group: {group}{condition}  ---")
print('organ number', len(df['Name'].unique()))
print('chr number', len(df['Chr_ID'].unique()))
print('Data shape:', df.shape)

size = 3

# ==============================
# PLOT 1: A2T
# ==============================
df_sorted_A2T = df.sort_values('A2T_P-A2T_C')
x = np.arange(len(df_sorted_A2T))

plt.figure(figsize=(20, 10))
plt.scatter(x, df_sorted_A2T['A2T_P'], color='g', s=size+1)
plt.scatter(x, df_sorted_A2T['A2T_C'], color='b', s=size+1)
plt.scatter(x, df_sorted_A2T['A2T_annot'], color='c', s=size)
plt.scatter(x, df_sorted_A2T['A2T_W'], color='black', s=size)
plt.title(f'A2T P and C {group}{condition}')
plt.ylabel('A2T')
plt.xlabel('Genomes')
plt.legend()
plt.show()

# Clear and close the figure completely
#plt.clf()
#plt.close()

# Print A2T Stats
print('A2T_P mean', df['A2T_P'].mean())
print('A2T_P std', df['A2T_P'].std())
print('A2T_C mean', df['A2T_C'].mean())
print('A2T_C std', df['A2T_C'].std())
print('A2T_annot mean', df['A2T_annot'].mean())
print('A2T_annot std', df['A2T_annot'].std())
print('A2T_W mean', df['A2T_W'].mean())
print('A2T_W std', df['A2T_W'].std())







#plot 2


df_sorted_G2C = df.sort_values('G2C_P-G2C_C')
x = np.arange(len(df_sorted_G2C))

plt.figure(figsize=(20, 10))
plt.scatter(x, df_sorted_G2C['G2C_P'], color='yellow', s=size+1)
plt.scatter(x, df_sorted_G2C['G2C_C'], color='r', s=size+1)
plt.scatter(x, df_sorted_G2C['G2C_annot'], color='blue', s=size)
plt.scatter(x, df_sorted_G2C['G2C_W'], color='black', s=size)
plt.title(f'G2C P and C {group}{condition}')
plt.ylabel('G2C')
plt.xlabel('Genomes')
plt.legend()
plt.show()


# Print G2C Stats
print('G2C_P mean', df['G2C_P'].mean())
print('G2C_P std', df['G2C_P'].std())
print('G2C_C mean', df['G2C_C'].mean())
print('G2C_C std', df['G2C_C'].std())
print('G2C_annot mean', df['G2C_annot'].mean())
print('G2C_annot std', df['G2C_annot'].std())
print('G2C_W mean', df['G2C_W'].mean())
print('G2C_W std', df['G2C_W'].std())





'''for group in groups:
    # Filter data for the current group
    df = my_datafram[my_datafram['Genome'] == group]

    print(f"--- Group: {group} ---")
    print('organ number', len(df['Name'].unique()))
    print('chr number', len(df['Chr_ID'].unique()))
    print('Data shape:', df.shape)

    # ==============================
    # PLOT 1: A2T
    # ==============================
    df_sorted_A2T = df.sort_values('A2T_P-A2T_C')
    x = np.arange(len(df_sorted_A2T))
    size = 0.5

    plt.figure(figsize=(20, 10))
    plt.scatter(x, df_sorted_A2T['A2T_P'], color='g', s=size+1)
    plt.scatter(x, df_sorted_A2T['A2T_C'], color='b', s=size+1)
    plt.scatter(x, df_sorted_A2T['A2T_annot'], color='c', s=size)
    plt.scatter(x, df_sorted_A2T['A2T_W'], color='black', s=size)
    plt.title(f'A2T P and C {group}')
    plt.ylabel('A2T')
    plt.xlabel('Genomes')
    plt.legend()
    plt.show()

    # Clear and close the figure completely
    #plt.clf()
    #plt.close()

    # Print A2T Stats
    print('A2T_P mean', df['A2T_P'].mean())
    print('A2T_P std', df['A2T_P'].std())
    print('A2T_C mean', df['A2T_C'].mean())
    print('A2T_C std', df['A2T_C'].std())
    print('A2T_annot mean', df['A2T_annot'].mean())
    print('A2T_annot std', df['A2T_annot'].std())
    print('A2T_W mean', df['A2T_W'].mean())
    print('A2T_W std', df['A2T_W'].std())

    # ==============================
    # PLOT 2: G2C
    # ==============================
    df_sorted_G2C = df.sort_values('G2C_P-G2C_C')
    x = np.arange(len(df_sorted_G2C))

    plt.figure(figsize=(20, 10))
    plt.scatter(x, df_sorted_G2C['G2C_P'], color='yellow', s=size+1)
    plt.scatter(x, df_sorted_G2C['G2C_C'], color='r', s=size+1)
    plt.scatter(x, df_sorted_G2C['G2C_annot'], color='purple', s=size)
    plt.scatter(x, df_sorted_G2C['G2C_W'], color='black', s=size)
    plt.title(f'G2C P and C {group}')
    plt.ylabel('G2C')
    plt.xlabel('Genomes')
    plt.legend()
    plt.show()

    # Clear and close the figure completely
    plt.clf()
    plt.close()

    # Print G2C Stats
    print('G2C_P mean', df['G2C_P'].mean())
    print('G2C_P std', df['G2C_P'].std())
    print('G2C_C mean', df['G2C_C'].mean())
    print('G2C_C std', df['G2C_C'].std())
    print('G2C_annot mean', df['G2C_annot'].mean())
    print('G2C_annot std', df['G2C_annot'].std())
    print('G2C_W mean', df['G2C_W'].mean())
    print('G2C_W std', df['G2C_W'].std())'''







