import os
import time
import psutil
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor
from Bio import SeqIO
from function_for_viruses import seq_Nucleotide_composition
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm  # Progress bar library

def process_directory(directory, path):
    if '.json' in directory:
        return None

    fasta_file = None
    gff_file = None
    files = os.listdir(path + '/' + directory)

    for file in files:
        if file.endswith('.fna'):
            fasta_file = file
        if file.endswith('.gff'):
            gff_file = file

    if not fasta_file:
        return None

    sequence = SeqIO.parse(path + '/' + directory + '/' + fasta_file, 'fasta')
    whole_genome = ''
    for record in sequence:
        whole_genome += str(record.seq)

    len_seqs = len(whole_genome)
    A, T, C, G = seq_Nucleotide_composition(whole_genome)

    return len_seqs, A, T, C, G, round(A / T, 2), round(G / C, 2)

def main():
    start = time.time()
    path = '/media/genomics/HDD/Genomics/Genomes/unzip_viruse_all_complete/ncbi_dataset/data'
    length_list = []
    A_list = []
    T_list = []
    C_list = []
    G_list = []
    A_T_list = []
    G_C_list = []

    list_direct = os.listdir(path)

    total_memory = psutil.virtual_memory().total
    memory_limit = total_memory   # 90% of total memory
    num_workers = max(1, int(os.cpu_count() ))  # 90% of available CPUs
    print(f"Using {num_workers} worker processes")

    progress = tqdm(total=len(list_direct), desc="Processing directories", unit="dir")  # Initialize progress bar

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = []
        for directory in list_direct:
            futures.append(executor.submit(process_directory, directory, path))

            if len(futures) >= num_workers:
                for future in concurrent.futures.as_completed(futures):
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

                    progress.update(1)  # Update progress bar

                futures = []  # Clear the list of futures to free memory

                # Monitor memory usage
                memory_info = psutil.virtual_memory()
                while memory_info.used >= memory_limit:
                    print(f"Memory usage {memory_info.percent}% reached. Pausing to free memory.")
                    time.sleep(1)
                    memory_info = psutil.virtual_memory()

        # Process remaining futures
        for future in concurrent.futures.as_completed(futures):
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

            progress.update(1)  # Update progress bar

    progress.close()  # Close the progress bar

    x = np.arange(1, len(length_list) + 1, 1)
    print(len(x), len(length_list))

    print(np.mean(length_list))
    plt.figure()
    plt.scatter(x, A_T_list, color='green', s=0.1)
    plt.scatter(x, G_C_list, color='red', s=0.1)
    plt.savefig('scatter_virus_parallel.jpg')
    plt.figure()
    plt.scatter(x,A_list)
    plt.scatter(x,T_list)
    plt.scatter(x,C_list)
    plt.scatter(x,G_list)
    plt.savefig('scatter_virus_A_T_C_G_parallel.jpg', )

    plt.figure()
    plt.scatter(x,length_list)
    plt.savefig('scatter_virus_length_parallel.jpg')

    end = time.time()
    print(f"Total time: {end - start}")

if __name__ == "__main__":
    main()
