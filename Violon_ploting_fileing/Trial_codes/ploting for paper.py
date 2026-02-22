import matplotlib.pyplot as plt
import pandas as pd

file_path = r'D:\PhD\chargaf_genome_test_project\fungi2_git\fungi\refrence_CDS_2025-01-07_23-21-41\results_2025-01-07_23-21-41.xlsx'
file = pd.read_excel(file_path)
import numpy as np

path_to_save = r'D:\PhD\chargaf_genome_test_project\fungi2_git\fungi\refrence_CDS_2025-01-07_23-21-41\\'


plt.scatter(file['GC_P'],file['P_A_usage_list'],color = 'green', s = 3 )
plt.scatter(file['GC_P'],file['P_T_usage_list'],color = 'blue',  s = 3 )
plt.scatter(file['GC_P'],file['P_C_usage_list'],color = 'yellow',s = 3 )
plt.scatter(file['GC_P'],file['P_G_usage_list'],color = 'red',   s = 3 )
plt.xlabel('GC content(%)')
plt.ylabel('nucleotide usage P %')
plt.title('Fungi\'s CDS nucleotide usage')
plt.savefig(fr'{path_to_save}\'Fungi CDS nucleotide usage P.svg',transparent = True,)
plt.show()

plt.scatter(file['GC_C'],file['C_A_usage_list'],color = 'green',s = 5 )
plt.scatter(file['GC_C'],file['C_T_usage_list'],color = 'blue',s = 5 )
plt.scatter(file['GC_C'],file['C_C_usage_list'],color = 'yellow',s = 5 )
plt.scatter(file['GC_C'],file['C_G_usage_list'],color = 'red',s = 5 )
plt.xlabel('GC content(%)')
plt.ylabel('nucleotide usage C %')
plt.title('Bacteria\'s CDS nucleotide usage')
plt.savefig(fr'{path_to_save}\'Bacterias CDS nucleotide usage C.svg',transparent = True,)
plt.show()

plt.scatter(file['Length(M)'],file['P2C'],color = 'black',s = 5 )
plt.xlabel('Length(M)')
plt.ylabel('P2C')
plt.title('P2C to length Bacteria')
plt.show()



plt.scatter(file['GC_P'],file['A2T_P'],color = 'green',s = 5 )
plt.scatter(file['GC_C'],file['A2T_C'],color = 'blue',s = 5 )
plt.scatter(file['GC_P'],file['C2G_P'],color = 'yellow',s = 5 )
plt.scatter(file['GC_C'],file['C2G_C'],color = 'red',s = 5 )
plt.xlabel('GC content(%)')
plt.ylabel('nucleotide usage C %')
plt.title('Bacteria\'s A2T and C2G in P and C')
plt.savefig(fr'{path_to_save}\'A2T and C2G in P and C.svg',transparent = True,)
plt.show()

plt.scatter(file['dif2sum'],file['GC_W'],color = 'green',s = 5 )

plt.xlabel('GC content(%)')
plt.ylabel('dif2sum')
plt.title('Fungis\'s dif2sum')
plt.savefig(fr'{path_to_save}\'dif2sum.svg',transparent = True,)
plt.show()


import pandas as pd
import numpy as np

# Example data structure, assuming 'data' is your DataFrame with columns 'genome_length' and 'P2C'
data = file

# Calculate the mean, standard deviation, correlation, and other statistics
mean_p2c = data['dif2sum'].mean()
std_p2c = data['dif2sum'].std()
correlation = data['Length(M)'].corr(data['dif2sum'])
min_p2c = data['dif2sum'].min()
max_p2c = data['dif2sum'].max()
range_p2c = max_p2c - min_p2c
median_p2c = data['dif2sum'].median()

# Print out the results
print(f"Mean P2C: {mean_p2c}")
print(f"Standard Deviation of P2C: {std_p2c}")
print(f"Correlation between genome length and P2C: {correlation}")
print(f"Range of P2C: {range_p2c}")
print(f"Median P2C: {median_p2c}")


