import time

start = time.time()
v_number = 10000
import  concurrent.futures
import time
def process_data(viruse_name):
    folder_path = '/media/genomics/HDD/Genomics/Genomes/unzip_viruse_all_complete/ncbi_dataset/data'

    from Bio import SeqIO
    from function_for_viruses import seq_Nucleotide_composition
    length_list = []
    A_list = []
    T_list = []
    C_list = []
    G_list = []
    A_T_list = []
    G_C_list = []

    files = os.listdir(folder_path+'/'+viruse_name)
    #print(files)
    for file in files:
        if file.endswith('.fna'):
            fasta_file = file
        if file.endswith('.gff'):
            gff_file = file
    sequence = SeqIO.parse(folder_path+'/'+viruse_name+'/'+fasta_file,'fasta')
    rec_counter = 0
    sequence_list = []
    for record in sequence:
        rec_counter +=1
        #print(f'recored{rec_counter}')
        whole_genome = str(record.seq)
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
        #print(['start',A_list,T_list,C_list,G_list,A_T_list,G_C_list,length_list,'end'])
    return [A_list,T_list,C_list,G_list,A_T_list,G_C_list,length_list]


import os
viruses_name = os.listdir('/media/genomics/HDD/Genomics/Genomes/unzip_viruse_all_complete/ncbi_dataset/data')
viruses_name = viruses_name[0:v_number]
#print(viruses_name)
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executer:
    results = list(executer.map(process_data,viruses_name))
end = time.time()
A = []
T = []
C = []
G = []
A_T = []
G_C = []
L = []

extended = []
for i in results  :
    A.extend(i[0])
    T.extend(i[1])
    C.extend(i[2])
    G.extend(i[3])
    A_T.extend(i[4])
    G_C.extend(i[5])
    L.extend(i[6])


import  numpy as np
import matplotlib.pyplot as plt
x = np.arange(1,len(A_T)+1,1)
#plt2.plot(A_T_list,'green')
#plt2.plot(G_C_list,'red')
#plt.savefig('P2virus3.svg',format ='svg')
print(np.mean(L))

plt.scatter(x,A_T,color = 'green',s = 0.1)
plt.scatter(x,G_C,color = 'red', s = 0.1)
plt.savefig('scater_virus_parrallel.svg',format ='svg')





print(len(A),len(T),len(C),len(G),len(A_T),len(G_C),len(L))

#end = time.time()
print(end - start)

print('time',end-start)