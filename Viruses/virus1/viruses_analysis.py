
path = r'/home/sukimi/unzip_viruse_all_complete/ncbi_dataset/data'
#path =  r'D:\PhD\genomes\unzip_viruse_all_complete\ncbi_dataset\data'
import os
from Bio import SeqIO
RNA_minis = ['Tospoviridae','Phenuiviridae','Jonviridae','Hantaviridae','Fimoviridae','Feraviridae','Tenuivirus','Filoviridae','Peribunyaviridae','Paramyxoviridae','Rhabdoviridae','Orthomyxoviridae','Varicosavirus','Ophiovirus','Arenaviridae']
RNA_plus = []
length_list = []

counter = 0
list_direct = os.listdir(path)
for i in list_direct:
    counter+=1
    if counter == 1000:
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
    #print('new viruse *******************************************************************************')
    for record in sequence:
        #print('record',record.id)

        sequence_segment = record.seq
        len_seqs = len(record.seq)
        #print(len_seqs)
        length_list.append(len_seqs)
        inf = record.description
        #print('ch',record.id)
        for name in RNA_minis:
            if name in inf :
                print('-RNA')
                gff = open(path+'/'+i+'/'+gff_file)
                for line in gff:
                    if '#' in line :
                        continue
                    line = line.strip()
                    line = line.split()
                    if line[0] == record.id:
                        print(line[0],record.id)
                    print(record.id)
                    print(line[6])
                    print(len(record.seq))
                    start = int(line[3])
                    end = int(line[4])
                    print(sequence_segment[start:end])
            #if name in

import matplotlib.pyplot as plt

import numpy as np
x = np.arange(1,len(length_list)+1,1)
print(len(x),len(length_list))
#plt.plot(length_list)
length_list.sort()
plt.scatter(x,length_list,color = 'green',s = 10)
#plt.show()
plt.savefig('virus1.svg',format = 'svg',dpi = 100,transparent=True)
plt.show()

print(np.mean(length_list))


