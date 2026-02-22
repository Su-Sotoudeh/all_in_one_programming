import time

start = time.time()
path = '/media/genomics/HDD/Genomics/Genomes/unzip_viruse_all_complete/ncbi_dataset/data'
import os
from Bio import SeqIO
from function_for_viruses import seq_Nucleotide_composition
length_list = []
A_list = []
T_list = []
C_list = []
G_list = []
A_T_list = []
G_C_list = []
counter = 0
list_direct = os.listdir(path)
for i in list_direct:
    counter+=1
    if counter == 10000:
        break
    if '.json' in i:
        continue
    files = os.listdir(path+'/'+i)
    #print(files)
    for file in files:
        if file.endswith('.fna'):
            fasta_file = file
        if file.endswith('.gff'):
            gff_file = file
    sequence = SeqIO.parse(path+'/'+i+'/'+fasta_file,'fasta')
    whole_genome = ''
    for record in sequence:
        whole_genome += str(record.seq)
    len_seqs = len(whole_genome)
    #print(len_seqs)
    A,T,C,G = seq_Nucleotide_composition(whole_genome)
    A_list.append(A)
    T_list.append(T)
    C_list.append(C)
    G_list.append(G)
    A_T_list.append(round(A/T,2))
    G_C_list.append(round(G/C,2))
    length_list.append(len_seqs)

import matplotlib.pyplot as plt

import numpy as np
x = np.arange(1,len(length_list)+1,1)
print(len(x),len(length_list))
#plt.plot(length_list)
#length_list.sort()
'''plt.scatter(length_list,A_list,color = 'green',s = 10)
plt.scatter(length_list,T_list,color = 'blue',s = 10)
plt.scatter(length_list,C_list,color = 'yellow',s = 10)
plt.scatter(length_list,G_list,color = 'red',s = 10)
plt.savefig('P2virus2.svg',format = 'svg')#transparent=True'''
import  numpy as np
import matplotlib.pyplot as plt2
x = np.arange(1,len(A_T_list)+1,1)
#plt2.plot(A_T_list,'green')
#plt2.plot(G_C_list,'red')
#plt.savefig('P2virus3.svg',format ='svg')
print(np.mean(length_list))

plt2.scatter(x,A_T_list,color = 'green',s = 0.1)
plt2.scatter(x,G_C_list,color = 'red', s = 0.1)
plt.savefig('scater_viruse_non_paralel.svg',format ='svg')


end = time.time()
print(end - start)


