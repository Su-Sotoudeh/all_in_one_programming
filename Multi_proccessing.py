'''import multiprocessing
import time
sequences = ['ATCGAGAGGGGAGGAAGTTGGGGG','AAAAAAAAAAAAAAAAAAAAAAA','TATAAATAGAGGGAGAGAGGAATAATATATATA','GGGGGGGG']
def worker(sequence,gene_number):
    A = sequence.count('A')
    T = sequence.count('T')
    C = sequence.count('C')
    G = sequence.count('G')
    #print(gene_number,'A:',A,'  T:',T,' C:',C,' G:',G)


if __name__ == "__main__":
    processes = []
    for index,seq in enumerate(sequences):
        print(index,seq)
        p = multiprocessing.Process(target=worker, args=(seq,index,))

        processes.append(p)
        print(processes)
        p.start()

    for p in processes:
        p.join()  # Wait for all processes to complete

    print("All processes completed")'''

'''import multiprocessing
import time

# Function to count nucleotides in a sequence
def count_nucleotides(sequence, name):
    counts = {'A': sequence.count('A'), 'T': sequence.count('T'),
              'C': sequence.count('C'), 'G': sequence.count('G')}
    return name, counts

# Sample sequence
sequence = "ATCG" * 1000  # 4 million bases

# Single-process version
def single_process():
    start = time.time()
    result = count_nucleotides(sequence, "seq1")
    print(f"Single process result: {result}")
    print(f"Time: {time.time() - start:.2f} seconds")

# Multiprocessing version
def multi_process():
    start = time.time()
    # Split sequence into 4 chunks
    chunk_size = len(sequence) // 4
    chunks = [sequence[i:i + chunk_size] for i in range(0, len(sequence), chunk_size)]
    #print(chunks)
    # Create a pool with 4 workers
    with multiprocessing.Pool(processes=4) as pool:

        results = pool.starmap(count_nucleotides, [(chunk, f"chunk{i}") for i, chunk in enumerate(chunks)])
        #print(results)
    # Combine results
    print(results)
    total_counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    for name, counts in results:
        for base in total_counts:
            total_counts[base] += counts[base]
    print(f"Multiprocessing result: {total_counts}")
    print(f"Time: {time.time() - start:.2f} seconds")

if __name__ == "__main__":
    print("Running single process...")
    single_process()
    print("\nRunning multiprocessing...")
    multi_process()'''


"""
simple_chr_parallel.py

Simple, minimal example: one worker per chromosome, counting A/T/C/G/N and writing
one CSV row per chromosome. Shows how multiprocessing works step-by-step.
"""









import os
import gzip
import time
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
import pandas as pd
from Bio import SeqIO

# --------- CONFIG ---------
DATA_DIR = r"D:\PhD\genomes\bacteria\ncbi_dataset\data\\"
import os
import gzip
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import SeqIO

# ------------ Worker function (runs in parallel) ------------
def process_chromosome(fasta_path, chr_id, organism_name):
    # open fasta file
    handle = gzip.open(fasta_path, "rt") if fasta_path.endswith(".gz") else open(fasta_path, "r")
    seq = None
    for rec in SeqIO.parse(handle, "fasta"):
        if rec.id == chr_id:
            seq = str(rec.seq)
            break
    handle.close()

    if seq is None:
        return None

    # simple counts
    A = seq.count("A") + seq.count("a")
    T = seq.count("T") + seq.count("t")
    C = seq.count("C") + seq.count("c")
    G = seq.count("G") + seq.count("g")

    # return one row (dictionary)
    return {
        "Name": organism_name,
        "Chr_ID": chr_id,
        "Length": len(seq),
        "A": A,
        "T": T,
        "C": C,
        "G": G,
    }

# ------------ Main driver (submits jobs & writes CSV) ------------
def main():
    data_dir = r"D:\PhD\genomes\bacteria\ncbi_dataset\data\\"   # <-- change this
    results_csv = "chr_counts.csv"

    tasks = []
    # collect all (organism, fasta, chr_id)
    for org in os.listdir(data_dir):
        org_path = os.path.join(data_dir, org)
        if not os.path.isdir(org_path):
            continue
        fasta_files = [f for f in os.listdir(org_path) if f.endswith(".fna") or f.endswith(".fna.gz")]
        if not fasta_files:
            continue
        fasta_path = os.path.join(org_path, fasta_files[0])
        # read chromosome IDs
        handle = gzip.open(fasta_path, "rt") if fasta_path.endswith(".gz") else open(fasta_path, "r")
        for rec in SeqIO.parse(handle, "fasta"):
            tasks.append((fasta_path, rec.id, org))
        handle.close()

    print(f"Found {len(tasks)} chromosomes to process.")

    results = []
    # run in parallel
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_chromosome, f, c, o) for f, c, o in tasks]
        for future in as_completed(futures):
            res = future.result()
            if res:
                results.append(res)
                print("Finished:", res["Name"], res["Chr_ID"])

    # save all results
    df = pd.DataFrame(results)
    df.to_csv(results_csv, index=False)
    print(f"Saved results to {results_csv}")

# ------------ Run ------------
if __name__ == "__main__":
    main()










