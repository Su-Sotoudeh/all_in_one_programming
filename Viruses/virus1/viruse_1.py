
path = r'D:\PhD\genomes\unzip_viruse_all_complete\ncbi_dataset\data'
import os
from Bio import SeqIO

length_list = []

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
    for record in sequence:
        len_seqs = len(record.seq)
        #print(len_seqs)
        length_list.append(len_seqs)

import matplotlib.pyplot as plt

import numpy as np
x = np.arange(1,len(length_list)+1,1)
print(len(x),len(length_list))
#plt.plot(length_list)
length_list.sort()
plt.scatter(x,length_list,color = 'green',s = 10)
plt.show()
plt.savefig('virus1.svg',format = 'svg',dpi = 100,transparent=True)

print(np.mean(length_list))


