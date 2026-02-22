import os
import time
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor
from Bio import SeqIO
from function_for_viruses import seq_Nucleotide_composition
import matplotlib.pyplot as plt
import numpy as np

def process_directory(i):
    if '.json' in i:
        return None

    fasta_file = None
    gff_file = None
    files = os.listdir(path + '/' + i)

    for file in files:
        if file.endswith('.fna'):
            fasta_file = file
        if file.endswith('.gff'):
            gff_file = file

    if not fasta_file:
        return None

    sequence = SeqIO.parse(path + '/' + i + '/' + fasta_file, 'fasta')
    whole_genome = ''
    for record in sequence:
        whole_genome += str(record.seq)

    len_seqs = len(whole_genome)
    A, T, C, G = seq_Nucleotide_composition(whole_genome)

    return len_seqs, A, T, C, G, round(A / T, 2), round(G / C, 2)

def main():
    start = time.time()
    global path
    path = '/media/genomics/HDD/Genomics/Genomes/unzip_viruse_all_complete/ncbi_dataset/data'
    length_list = []
    A_list = []
    T_list = []
    C_list = []
    G_list = []
    A_T_list = []
    G_C_list = []

    list_direct = os.listdir(path)

    num_workers = os.cpu_count()  # Use the number of logical CPUs
    print(f"Using {num_workers} worker processes")

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(process_directory, i) for i in list_direct[:]]

        for future in futures:
            result = future.result()
            if result:
                len_seqs, A, T, C, G, A_T, G_C = result
                length_list.append(len_seqs)
                A_list.append(A)
                T_list.append(T)
                C_list.append(C)
                G_list.append(G)
                A_T_list.append(A_T)
                G_C_list.append(G_C)

    x = np.arange(1, len(length_list) + 1, 1)
    print(len(x), len(length_list))

    print(np.mean(length_list))

    plt.scatter(x, A_T_list, color='green', s=0.1)
    plt.scatter(x, G_C_list, color='red', s=0.1)
    plt.savefig('scater_viruse_non_paralel.svg', format='svg')

    end = time.time()
    print(end - start)

if __name__ == "__main__":
    main()
