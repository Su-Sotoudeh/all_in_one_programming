import os
import time
import psutil
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor
from Bio import SeqIO
#from function_for_viruses import seq_Nucleotide_composition

from function_for_viruses import calculat_portion_of_elements_overlap_sharing
from function_for_viruses import overlap_nucleotide_composition

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from tqdm import tqdm  # Progress bar library

def process_directory(directory, path):
    if '.json' in directory:
        return None

    fasta_file = None
    gff_file = None
    files = os.listdir(os.path.join(path, directory))
    print(files)
    for file in files:
        if file.endswith('.fna'):
            fasta_file = file
        if file.endswith('.gff'):
            gff_file = file
            break
    if not fasta_file:
        return None

    sequence = SeqIO.parse(os.path.join(path, directory, fasta_file), 'fasta')
    whole_genome = ''

    for record in sequence:

        whole_genome = str(record.seq)

        len_seqs = len(whole_genome)


        labled_sequence = calculat_portion_of_elements_overlap_sharing(whole_genome,os.path.join(path, directory, gff_file),record.id,['CDS'])
        if labled_sequence == 'map error' :
            continue

        result = overlap_nucleotide_composition(whole_genome)
        if result == 'unsequenced_region > 20':
            continue
        A,T,C,G,L = result


        return len_seqs, A, T, C, G, round(A / T, 2), round(G / C, 2)

def main():
    start = time.time()
    path = r'D:\PhD\genomes\unzip_viruse_all_complete\ncbi_dataset\data'
    length_list = []
    A_list = []
    T_list = []
    C_list = []
    G_list = []
    A_T_list = []
    G_C_list = []

    list_direct = os.listdir(path)
    #list_direct = list_direct[:1000]

    total_memory = psutil.virtual_memory().total
    memory_limit = total_memory   # 90% of total memory
    num_workers = max(1, int(os.cpu_count()))  # 90% of available CPUs
    print(f"Using {num_workers} worker processes")

    progress = tqdm(total=len(list_direct), desc="Processing directories", unit="dir")  # Initialize progress bar

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = []
        counter = 0
        for directory in list_direct:
            counter += 1
            if counter == 100:
                break
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

    sns.set(style="whitegrid")  # Set the Seaborn style
    x= length_list
    plt.figure()
    sns.scatterplot(x=x, y=A_T_list, color='green', s=10, label='A-T')
    sns.scatterplot(x=x, y=G_C_list, color='red', s=10, label='G-C')
    plt.legend()
    plt.xlabel('Index')
    plt.ylabel('Ratio')
    plt.title('A-T and G-C Ratios')
    plt.savefig('scatter_virus_parallel.jpg')
    #plt.show()

    plt.figure()
    sns.scatterplot(x=x, y=A_list, color='blue', s=10, label='A')
    sns.scatterplot(x=x, y=T_list, color='orange', s=10, label='T')
    sns.scatterplot(x=x, y=C_list, color='purple', s=10, label='C')
    sns.scatterplot(x=x, y=G_list, color='brown', s=10, label='G')
    plt.legend()
    plt.xlabel('Index')
    plt.ylabel('Nucleotide Count')
    plt.title('Nucleotide Counts')
    plt.savefig('scatter_virus_A_T_C_G_parallel.jpg')
    plt.show()

    plt.figure()
    sns.scatterplot(x=x, y=length_list, color='cyan', s=10)
    plt.xlabel('Index')
    plt.ylabel('Genome Length')
    plt.title('Genome Lengths')
    plt.savefig('scatter_virus_length_parallel.jpg')
    plt.show()

    end = time.time()
    print(f"Total time: {end - start}")

if __name__ == "__main__":
    main()