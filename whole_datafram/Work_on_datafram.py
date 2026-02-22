import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
df = pd.read_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\whole_datafram\Whole_datafram_2026.csv',low_memory=False)
df = pd.read_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\Viruses\Viruses_refrence_Int_UTR_CDS_2026-01-26_16-51-15\Viruses_2026-01-26_16-51-15.csv',low_memory=False)
viruse_type = 'ssRNA(-)'
df = df[df['Genome'] == viruse_type]
df = df[df['A2T_P'] != 0]
df = df[df['G2C_P'] != 0]

groups = df['Group'].unique()

#df =df[(df['Chr_ID'].str.startswith('NC'))&(df['Wh_ano']>10) &(df['A2T_W']>0.98)&(df['A2T_W']<1.02)]
print(df.shape)
df = df.sort_values('A2T_P-A2T_C')

A2T_P = df['A2T_P']
A2T_C = df['A2T_C']
A2T_annotation = df['A2T_annot']
A2T_W = df['A2T_W']
x = np.arange(len(A2T_P))
size = 0.5
# Example: Width 20 inches, Height 10 inches
fig, axes = plt.subplots( figsize=(20, 10))
plt.scatter(x,A2T_P,color='g',s = size)
plt.scatter(x,A2T_C,color='r',s = size)
#plt.scatter(x,A2T_annotation,color='b',s=size)
#plt.scatter(x,A2T_annotation,color='y',s=size)
#plt.scatter(x,A2T_W,color='black',s=size)
plt.title(f'A2T P and C {viruse_type}')
plt.ylabel('A2T')
plt.xlabel('Genomes')
plt.legend()
plt.show()
print('A2T_P mean',df['A2T_P'].mean())
print('A2T_P std',df['A2T_P'].std())
print('A2T_C mean',df['A2T_C'].mean())
print('A2T_C std',df['A2T_C'].std())
print('A2T_annot mean',df['A2T_annot'].mean())
print('A2T_annot std',df['A2T_annot'].std())
print('A2T_W mean',df['A2T_W'].mean())
print('A2T_W std',df['A2T_W'].std())



df = df.sort_values('G2C_P-G2C_C')

A2T_P = df['G2C_P']
A2T_C = df['G2C_C']
A2T_annotation = df['G2C_annot']
A2T_W = df['G2C_W']
x = np.arange(len(A2T_P))
size = 0.5
fig, axes = plt.subplots( figsize=(20, 10))
plt.scatter(x,A2T_P,color='g',s = size)
plt.scatter(x,A2T_C,color='r',s = size)
#plt.scatter(x,A2T_annotation,color='b',s=size)
#plt.scatter(x,A2T_annotation,color='y',s=size)
#plt.scatter(x,A2T_W,color='black',s=size)

plt.title(f'G2C P and C {viruse_type}')
plt.ylabel('G2C')
plt.xlabel('Genomes')
plt.legend()
plt.show()

print('G2C_P mean',df['G2C_P'].mean())
print('G2C_P std',df['G2C_P'].std())
print('G2C_C mean',df['G2C_C'].mean())
print('G2C_C std',df['G2C_C'].std())
print('G2C_annot mean',df['G2C_annot'].mean())
print('G2C_annot std',df['G2C_annot'].std())
print('G2C_W mean',df['G2C_W'].mean())
print('G2C_W std',df['G2C_W'].std())


