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
DATA_DIR =   r"/media/su/Bioinf/Genomes_2025_3/Genomes_2025_3_7/Fungi_annotated_compelet_chromosome_genome_gff_cds_2025/ncbi_dataset/data/" # pc 5118

# <-- set this (organism subfolders expected)
RESULTS_CSV = f"chr_nt_counts_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
NUM_WORKERS = None   # None -> use all CPU cores
# If you want quick tests, set TEST_LIMIT to a small number (e.g. 5). Set to None to run all.
TEST_LIMIT = None

# --------- Worker: process one chromosome ---------
def process_chromosome(fasta_path, chr_id, organism_name=None):
    """
    Worker function executed in a separate process.
    Inputs:
      - fasta_path: path to .fna(.gz) file
      - chr_id: sequence id to process (string)
      - organism_name: optional, used for output labeling
    Returns:
      - dict with result fields (will be written to CSV by main process)
    """
    proc_name = multiprocessing.current_process().name
    # Log start (stdout from worker; may interleave with main prints)
    print(f"{proc_name}: starting {organism_name or ''}:{chr_id}")

    # open fasta (support gzipped fasta)
    handle = gzip.open(fasta_path, "rt") if fasta_path.endswith(".gz") else open(fasta_path, "r")
    seq_record = None
    for rec in SeqIO.parse(handle, "fasta"):
        if rec.id == chr_id:
            seq_record = rec
            break
    handle.close()

    if seq_record is None:
        # Could not find the requested chromosome id
        return {"Name": organism_name or "", "Chr_ID": chr_id, "error": "not_found"}

    seq = str(seq_record.seq)
    # simple counts (case-insensitive)
    A = seq.count("A") + seq.count("a")
    T = seq.count("T") + seq.count("t")
    C = seq.count("C") + seq.count("c")
    G = seq.count("G") + seq.count("g")
    N = seq.count("N") + seq.count("n")
    L = len(seq) - N if (len(seq) - N) > 0 else len(seq)

    GC_pct = round((G + C) / (A + T + C + G) * 100, 4) if (A + T + C + G) != 0 else 0
    A2T = round(A / T, 4) if T != 0 else 0
    G2C = round(G / C, 4) if C != 0 else 0

    result = {
        "Name": organism_name or "",
        "Chr_ID": chr_id,
        "Length": L,
        "A": A,
        "T": T,
        "C": C,
        "G": G,
        "N": N,
        "GC_pct": GC_pct,
        "A2T": A2T,
        "G2C": G2C,
    }

    print(f"{proc_name}: finished {organism_name or ''}:{chr_id}")
    return result

# --------- Helper: enumerate chromosome tasks ---------
def list_chromosome_tasks(data_dir, test_limit=None):
    """
    Walk DATA_DIR, find organism subfolders, find fasta (.fna or .fna.gz).
    Yield tuples (organism_name, fasta_path, chr_id).
    """
    tasks = []
    for org in sorted(os.listdir(data_dir)):
        org_path = os.path.join(data_dir, org)
        if not os.path.isdir(org_path):
            continue
        # find fasta file(s) in organism folder
        fna_files = [f for f in os.listdir(org_path) if f.endswith(".fna") or f.endswith(".fna.gz")]
        if not fna_files:
            continue
        # choose the first fasta (adjust logic as you need)
        fasta_path = os.path.join(org_path, sorted(fna_files)[0])
        # read headers and create tasks (we only read headers here; this will also read sequences briefly)
        handle = gzip.open(fasta_path, "rt") if fasta_path.endswith(".gz") else open(fasta_path, "r")
        count = 0
        for rec in SeqIO.parse(handle, "fasta"):
            tasks.append((org, fasta_path, rec.id))
            count += 1
            if test_limit and count >= test_limit:
                break
        handle.close()
        # optional global test_limit controlling total tasks
        if test_limit and len(tasks) >= test_limit:
            break
    return tasks

# --------- Main controller (submits chromosome jobs) ---------
def main(data_dir, results_csv, num_workers=None, test_limit=None):
    if num_workers is None:
        num_workers = multiprocessing.cpu_count()
    num_workers = int(num_workers)

    # build tasks
    tasks = list_chromosome_tasks(data_dir, test_limit=test_limit)
    print(f"Found {len(tasks)} chromosome tasks. Using {num_workers} workers.")

    # resume logic: read already processed Chr_ID
    processed = set()
    if os.path.exists(results_csv):
        try:
            df_done = pd.read_csv(results_csv, usecols=["Chr_ID"])
            processed = set(df_done["Chr_ID"].astype(str).tolist())
            print(f"Resuming: {len(processed)} chromosomes already in {results_csv}")
        except Exception:
            processed = set()

    # only submit tasks not already processed
    tasks_to_run = [(org, fasta, ch) for (org, fasta, ch) in tasks if str(ch) not in processed]
    print(f"{len(tasks_to_run)} tasks to run (skipping already processed).")

    # if nothing to do:
    if not tasks_to_run:
        print("No new tasks. Exiting.")
        return

    # Main will be the only writer (safe). Submit jobs to worker processes.
    header_written = os.path.exists(results_csv)  # write header only if file absent
    submitted = 0
    start_time = time.time()

    with ProcessPoolExecutor(max_workers=num_workers) as exe:
        future_to_task = {}
        for org, fasta_path, ch_id in tasks_to_run:
            future = exe.submit(process_chromosome, fasta_path, ch_id, org)
            future_to_task[future] = (org, ch_id)
            submitted += 1

        # collect results as they finish
        for future in as_completed(future_to_task):
            org, ch_id = future_to_task[future]
            try:
                res = future.result()
            except Exception as exc:
                print(f"Task {org}:{ch_id} generated an exception: {exc}")
                continue

            # append to CSV (main process only)
            if res is None:
                continue
            df_row = pd.DataFrame([res])
            df_row.to_csv(results_csv, mode="a", header=not header_written, index=False)
            header_written = True
            print(f"Saved row for {org}:{ch_id} -> {results_csv}")

    elapsed = time.time() - start_time
    print(f"All done. Submitted {submitted} tasks. Elapsed {elapsed:.1f}s.")


# --------- run (protect for Windows) ---------
if __name__ == "__main__":
    # quick test: set TEST_LIMIT to a small number above or set to None to process everything
    main(DATA_DIR, RESULTS_CSV, num_workers=NUM_WORKERS, test_limit=TEST_LIMIT)













