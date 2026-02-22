import pandas as pd
import os
the_matric = []
for row in range(18):
    columns = [[],[],[]]
    the_matric.append(columns)

not_corect_files = ['Mus']
folders = os.listdir(r"C:\Genomics\Chargaff\coding\all_in_one_programming\CSV_output")
folder_counter =0
for folder in folders:

    if folder in not_corect_files:
        continue
    folder_counter += 1
    #if folder_counter == 2:
    #    break

    files = os.listdir(r"C:\Genomics\Chargaff\coding\all_in_one_programming\CSV_output" + '\\' + folder)
    chr_counter = 0
    for chr in files:
        chr_counter += 1
        #if chr_counter == 3 :
        #    break

        df = pd.read_csv(r"C:\Genomics\Chargaff\coding\all_in_one_programming\CSV_output" + '\\' + folder + '\\' + chr,nrows = 18)
        df.columns.values[0] = 'Element'
        nonsequenced = df.loc[16,'A']
        df.loc[16, 'A'] = nonsequenced/4
        df.loc[16, 'T'] = nonsequenced/4
        df.loc[16, 'C'] = nonsequenced/4
        df.loc[16, 'G'] = nonsequenced/4
        #this is for each element
        dic_new_df = {}
        names = []
        if int(df.iloc[17,5])/1000000 < 10:
            continue
        for m_row in range(18):
            names.append(df.iloc[m_row,0][:-1])
            try:
                A2T = round(float(df.iloc[m_row,1])/float(df.iloc[m_row,2]),3)
                #if A2T > 1.05 or A2T < 0.95:
                #print(folder,chr,df.iloc[m_row,0][:-1],A2T,int(df.iloc[17,5])/1000000)

                the_matric[m_row][0].append(A2T)
            except ZeroDivisionError:
                the_matric[m_row][0].append(1)
                #print(df.iloc[i,1]+df.iloc[7+i,1]),float(df.iloc[i,2]+df.iloc[7+i,2])
                #print(df.to_string())
            try:
                G2C = round(float(df.iloc[m_row,4])/float(df.iloc[m_row,3]),3)
                the_matric[m_row][1].append(G2C)
            except ZeroDivisionError:
                the_matric[m_row][1].append(1)
            the_matric[m_row][2].append(int(df.iloc[17,5]))



import matplotlib.pyplot as plt
import numpy as np

# Example 11x3x3 matrix (you can use your real data)
the_matrix = the_matric

# Row labels
row_labels = [
    "CDS+", "UTR+", "Intron", "LINE+", "SINE+", "LTR+",
    "DNA+","CDS-", "UTR-", "Intron-", "LINE-", "SINE-", "LTR-",
    "DNA-", "Shared", "Unannotated", "unsequenced", "whole"
]

# Calculate mean and std for first two columns
mean_col1 = [np.mean(row[0]) for row in the_matrix]
std_col1  = [np.std(row[0]) for row in the_matrix]

mean_col2 = [np.mean(row[1]) for row in the_matrix]
std_col2  = [np.std(row[1]) for row in the_matrix]

# X-axis positions
x = np.arange(len(the_matrix))  # 0 to 10
width = 0.35

# Plot
fig, ax = plt.subplots(figsize=(12, 6))

bars1 = ax.bar(x - width/2, mean_col1, width, yerr=std_col1, capsize=5, label='A2T')
bars2 = ax.bar(x + width/2, mean_col2, width, yerr=std_col2, capsize=5, label='G2C')

# X-axis labels
ax.set_xticks(x)
ax.set_xticklabels(row_labels, rotation=45, ha='right')

# Titles and labels
ax.set_xlabel('Element')
ax.set_ylabel('A/T G/C')
ax.set_title('all mammals')
ax.legend()
plt.ylim(0.4,1.7)
plt.tight_layout()
plt.savefig('Plot_E2_all_mammals.png')
plt.show()

