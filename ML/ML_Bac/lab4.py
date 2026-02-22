import gzip

import matplotlib
from Bio import SeqIO

from  func_for_bac import informational_composition_of_organism
from func_for_bac import seq_Nucleotide_composition
from func_for_bac import my_function_for_sorting_two_list
from func_for_bac import equality_usage_of_strand
import matplotlib.pyplot as plt


A_compusition_list = []
T_compusition_list = []
C_compusition_list = []
G_compusition_list = []
A_T_compusition_list  = [] # A minis T
G_C_compusition_list = [] # G minis C

GC_countent_list = []
GC_countent_list_w = []
len_seq_list = []
portion_of_maping_list = []
plus_and_minis_strand_portion = []



element = 'CDS'
A = r'D:\PhD\genomes\bacteria\ncbi_dataset_bac _refrence_annot.zip'
last_shered_foldur = r'D:\PhD\genomes\bacteria\ncbi_dataset\data'
import os
counter = 0
for file_name in os.listdir(last_shered_foldur):
    counter+=1
    '''if counter == 10:
        break'''
    print(counter)
    bac_foldur = fr'{last_shered_foldur}\{file_name}'
    file_list = []
    for file in os.listdir(bac_foldur):
        files = fr'{bac_foldur}\{file}'
        file_list.append(files)
    #print(file_list)
    if len(file_list) < 2:
        continue
    if file_list[0].endswith('.fna'):
        sequence_path = file_list[0]

    if file_list[1].endswith('.gff'):
        gff_file_path = file_list[1]

    for recoed in SeqIO.parse(sequence_path,'fasta') :
        if  not 'plasmid' in recoed.description :
            #print(recoed.description)
            #print('main genoeme')
            #print(recoed.id)
            seq = recoed.seq
            #print(len(seq))
            break
    #    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



    Aw,Tw,Cw,Gw = seq_Nucleotide_composition(seq)
    GC_countent_list_w.append(Gw+Cw)


    Plus, Comp, Both = informational_composition_of_organism(seq, gff_file_path, [element], 50000)

    A, T, C, G = seq_Nucleotide_composition(Both[2])

    A_compusition_list.append(A)
    T_compusition_list.append(T)
    C_compusition_list.append(C)
    G_compusition_list.append(G)
    A_T_compusition_list.append(A-T)
    G_C_compusition_list.append(G-C)
    # not true, GC_countent_list.append(G + C)
    GC_countent_list.append(G + C)
    len_seq_list.append(len(seq))
    portion_of_maping_list.append(len(Both[2])/len(seq)*100)
    '''if len(Both[2])> len(seq):
        print(Both[2])
        print('Hi')'''
    #2024/10/18
    plus,comp = equality_usage_of_strand(seq,gff_file_path)
    plus_and_minis_strand_portion.append(plus/comp)

import numpy as np
A = np.arange(0,len(plus_and_minis_strand_portion),1)
print(len(A),len(plus_and_minis_strand_portion))
plt.scatter(A,plus_and_minis_strand_portion)
plt.savefig('strand usage_equality.svg',transparent=True)
plt.show()


#new code 2024/06/22
A_T_sort,G_C_sort = my_function_for_sorting_two_list(A_T_compusition_list,G_C_compusition_list)
plt.scatter(GC_countent_list_w,A_compusition_list ,color='green',s = 10)
plt.scatter(GC_countent_list_w,T_compusition_list ,color='blue',s  = 10)
plt.scatter(GC_countent_list_w,C_compusition_list ,color='yellow',s= 10)
plt.scatter(GC_countent_list_w,G_compusition_list ,color='red', s  = 10)
plt.title('Bacteria ' + element)

plt.xlabel('GC %')
plt.ylabel('A,T,C,G %')
plt.savefig('1.svg',transparent=True)
plt.show()

plt.scatter(GC_countent_list_w,A_T_compusition_list, color = 'green')
plt.scatter(GC_countent_list_w,G_C_compusition_list, color = 'red')
plt.title('Bacteria')
plt.xlabel('GC%')
plt.ylabel('A-T and G-C')
plt.savefig('2.svg',transparent=True)
plt.show()

plt.plot(A_T_sort, color = 'green')
plt.plot(G_C_sort, color = 'red')
plt.xlabel('max to min diffrence')
plt.ylabel('A-T and G-C')
plt.savefig('3.svg',transparent=True)
plt.show()

plt.scatter(len_seq_list,GC_countent_list_w, color = 'black')
plt.title('Bacteria')
plt.xlabel('len_sequence')
plt.ylabel('GC %')
plt.savefig('4.svg',transparent=True)
plt.show()

plt.scatter(portion_of_maping_list,GC_countent_list_w)
plt.title('Bacteria')
plt.xlabel('annotated %')
plt.ylabel('GC %')
plt.savefig('5.svg',transparent=True)
plt.show()