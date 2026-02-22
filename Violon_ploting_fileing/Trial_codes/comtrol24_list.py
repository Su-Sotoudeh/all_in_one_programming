import os
import gzip
import pandas as pd
import numpy as np
import time
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed, wait, FIRST_COMPLETED
import multiprocessing
import concurrent.futures  # kept for clarity if you prefer qualified calls
from Bio import SeqIO
import psutil
from typing import List, Tuple
import heapq
import sys

# Verify Python version for concurrent.futures
if sys.version_info < (3, 3):
    raise ImportError("Python 3.3 or higher is required for concurrent.futures")

# ----------------- Constants for Memory Management -----------------
RAM_LIMIT_GB = 10
USABLE_FRACTION = 0.7  # Use 70% of RAM to leave buffer
BYTES_PER_BASE = 1.5  # Estimated bytes per base (bytearray + processing overhead)
PROCESS_OVERHEAD = 100 * 1024 * 1024  # 100MB per process for GFF/pandas
BUFFER_FACTOR = 1.2  # Require 20% extra free RAM
MAX_CORES = max(1, multiprocessing.cpu_count()-1 )  # leave one core free
BATCH_SIZE = 10  # Number of chromosomes to estimate per batch

# ----------------- Your Original Functions -----------------
def nucleotide_analysis_for_exon(sequence_plus):
    A_entron_P = sequence_plus.count('Q')
    T_entron_P = sequence_plus.count('W')
    C_entron_P = sequence_plus.count('E')
    G_entron_P = sequence_plus.count('R')
    sum_entron_P = A_entron_P + T_entron_P + C_entron_P + G_entron_P
    A_entron_C = sequence_plus.count('q')
    T_entron_C = sequence_plus.count('w')
    C_entron_C = sequence_plus.count('e')
    G_entron_C = sequence_plus.count('r')
    sum_entron_C = A_entron_C + T_entron_C + C_entron_C + G_entron_C
    A_exon_P = sequence_plus.count('S')
    T_exon_P = sequence_plus.count('D')
    C_exon_P = sequence_plus.count('F')
    G_exon_P = sequence_plus.count('H')
    sum_exon_P = A_exon_P + T_exon_P + C_exon_P + G_exon_P
    A_exon_C = sequence_plus.count('s')
    T_exon_C = sequence_plus.count('d')
    C_exon_C = sequence_plus.count('f')
    G_exon_C = sequence_plus.count('h')
    sum_exon_C = A_exon_C + T_exon_C + C_exon_C + G_exon_C

    A_CDS_P = sequence_plus.count('Z')
    T_CDS_P = sequence_plus.count('X')
    C_CDS_P = sequence_plus.count('V')
    G_CDS_P = sequence_plus.count('B')
    sum_CDS_P = A_CDS_P + T_CDS_P + C_CDS_P + G_CDS_P

    A_CDS_C = sequence_plus.count('z')
    T_CDS_C = sequence_plus.count('x')
    C_CDS_C = sequence_plus.count('v')
    G_CDS_C = sequence_plus.count('b')
    sum_CDS_C = A_CDS_C + T_CDS_C + C_CDS_C + G_CDS_C

    A_shared = sequence_plus.count('A')
    T_shared = sequence_plus.count('T')
    C_shared = sequence_plus.count('C')
    G_shared = sequence_plus.count('G')
    A_unannotated = sequence_plus.count('a')
    T_unannotated = sequence_plus.count('t')
    C_unannotated = sequence_plus.count('c')
    G_unannotated = sequence_plus.count('g')
    sum_unannotated = A_unannotated + T_unannotated + C_unannotated + G_unannotated

    A = [A_entron_P, A_entron_C, A_exon_P, A_exon_C, A_CDS_P, A_CDS_C, A_shared, A_unannotated]
    T = [T_entron_P, T_entron_C, T_exon_P, T_exon_C, T_CDS_P, T_CDS_C, T_shared, T_unannotated]
    C = [C_entron_P, C_entron_C, C_exon_P, C_exon_C, C_CDS_P, C_CDS_C, C_shared, C_unannotated]
    G = [G_entron_P, G_entron_C, G_exon_P, G_exon_C, G_CDS_P, G_CDS_C, G_shared, G_unannotated]
    dic = {'A': A, 'T': T, 'C': C, 'G': G}
    df = pd.DataFrame(dic)
    df.index = ['+entron', '-entron', '+UTR', '-UTR', '+CDS', '-CDS', 'Shared', 'Unannotated']
    return df

def exon_and_intron_masking(sequence, df_all, chr_id):
    ch_id_split = chr_id.split('.')[0]
    sequence = sequence.lower()
    sequence_plus = bytearray(sequence.encode('ascii'))
    Inreon_P_counter = 0
    Inreon_C_counter = 0
    UTR_P_counter = 0
    UTR_C_counter = 0
    CDS_P_counter = 0
    CDS_C_counter = 0

    df_all.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    filtered_df_all = df_all[df_all['seqid'].str.startswith(ch_id_split)]

    # gene (intron) masking
    filtered_df = filtered_df_all[filtered_df_all['type'].str.startswith('gene')]
    for annotation in filtered_df.to_numpy():
        start = int(annotation[3])
        end = int(annotation[4])
        if end - start == 0:
            continue
        if annotation[6] == '+':
            Inreon_P_counter += 1
            sl = sequence_plus[(start - 1):end]
            for i in range(len(sl)):
                ch = sl[i]
                if ch == ord('a'): sl[i] = ord('Q')
                elif ch == ord('t'): sl[i] = ord('W')
                elif ch == ord('c'): sl[i] = ord('E')
                elif ch == ord('g'): sl[i] = ord('R')
                elif ch == ord('q'): sl[i] = ord('A')
                elif ch == ord('w'): sl[i] = ord('T')
                elif ch == ord('e'): sl[i] = ord('C')
                elif ch == ord('r'): sl[i] = ord('G')
            sequence_plus[(start - 1):end] = sl
        if annotation[6] == '-':
            Inreon_C_counter += 1
            sl = sequence_plus[(start - 1):end]
            for i in range(len(sl)):
                ch = sl[i]
                if ch == ord('a'): sl[i] = ord('q')
                elif ch == ord('t'): sl[i] = ord('w')
                elif ch == ord('c'): sl[i] = ord('e')
                elif ch == ord('g'): sl[i] = ord('r')
                elif ch == ord('Q'): sl[i] = ord('A')
                elif ch == ord('W'): sl[i] = ord('T')
                elif ch == ord('E'): sl[i] = ord('C')
                elif ch == ord('R'): sl[i] = ord('G')
            sequence_plus[(start - 1):end] = sl

    # exon masking (uses UTR counters in this code)
    filtered_df = filtered_df_all[filtered_df_all['type'].str.startswith('exon')]
    for annotation in filtered_df.to_numpy():
        start = int(annotation[3])
        end = int(annotation[4])
        if end - start == 0:
            continue
        if annotation[6] == '+':
            UTR_P_counter += 1
            sl = sequence_plus[(start - 1):end]
            for i in range(len(sl)):
                ch = sl[i]
                if ch == ord('Q'): sl[i] = ord('S')
                elif ch == ord('W'): sl[i] = ord('D')
                elif ch == ord('E'): sl[i] = ord('F')
                elif ch == ord('R'): sl[i] = ord('H')
                elif ch == ord('q'): sl[i] = ord('A')
                elif ch == ord('w'): sl[i] = ord('T')
                elif ch == ord('e'): sl[i] = ord('C')
                elif ch == ord('r'): sl[i] = ord('G')
            sequence_plus[(start - 1):end] = sl
        if annotation[6] == '-':
            UTR_C_counter += 1
            sl = sequence_plus[(start - 1):end]
            for i in range(len(sl)):
                ch = sl[i]
                if ch == ord('q'): sl[i] = ord('s')
                elif ch == ord('w'): sl[i] = ord('d')
                elif ch == ord('e'): sl[i] = ord('f')
                elif ch == ord('r'): sl[i] = ord('h')
                elif ch == ord('Q'): sl[i] = ord('A')
                elif ch == ord('W'): sl[i] = ord('T')
                elif ch == ord('E'): sl[i] = ord('C')
                elif ch == ord('R'): sl[i] = ord('G')
            sequence_plus[(start - 1):end] = sl

    # CDS masking
    filtered_df = filtered_df_all[filtered_df_all['type'].str.startswith('CDS')]
    for annotation in filtered_df.to_numpy():
        start = int(annotation[3])
        end = int(annotation[4])
        if end - start == 0:
            continue
        if annotation[6] == '+':
            CDS_P_counter += 1
            sl = sequence_plus[(start - 1):end]
            for i in range(len(sl)):
                ch = sl[i]
                if ch == ord('S'): sl[i] = ord('Z')
                elif ch == ord('D'): sl[i] = ord('X')
                elif ch == ord('F'): sl[i] = ord('V')
                elif ch == ord('H'): sl[i] = ord('B')
                elif ch == ord('z'): sl[i] = ord('A')
                elif ch == ord('x'): sl[i] = ord('T')
                elif ch == ord('v'): sl[i] = ord('C')
                elif ch == ord('b'): sl[i] = ord('G')
            sequence_plus[(start - 1):end] = sl
        if annotation[6] == '-':
            CDS_C_counter += 1
            sl = sequence_plus[(start - 1):end]
            for i in range(len(sl)):
                ch = sl[i]
                if ch == ord('s'): sl[i] = ord('z')
                elif ch == ord('d'): sl[i] = ord('x')
                elif ch == ord('f'): sl[i] = ord('v')
                elif ch == ord('h'): sl[i] = ord('b')
                elif ch == ord('Z'): sl[i] = ord('A')
                elif ch == ord('X'): sl[i] = ord('T')
                elif ch == ord('V'): sl[i] = ord('C')
                elif ch == ord('B'): sl[i] = ord('G')
            sequence_plus[(start - 1):end] = sl

    sequence_plus = sequence_plus.decode('ascii')
    # Return order matches what process_one_chromosome expects in your code:
    return sequence_plus, Inreon_P_counter, Inreon_C_counter, UTR_P_counter, UTR_C_counter, CDS_P_counter, CDS_C_counter

# ----------------- Chromosome Processing Function -----------------
def process_one_chromosome(final_number, chr_counter, organism, chr_id, seq, df_all, all_genomes_folder_path):
    """
    Process a single chromosome (seq_record) and return its DataFrame row.
    Args:
        final_number: Total number of organisms (for printing).
        chr_counter: Global chromosome counter (for printing).
        organism: Organism name.
        chr_id: Chromosome ID.
        seq: Chromosome sequence (str).
        df_all: GFF DataFrame for the organism.
        all_genomes_folder_path: Base path (for context).
    Returns: pd.DataFrame with one row or None if skipped.
    """
    try:
        # Apply your original filters
        if len(seq) > 300000000:
            print(f"Skipping chromosome {chr_id} (organism: {organism}) due to size > 300M bases")
            return None
        l = len(seq) / 1000000
        if l < 0.05:
            print(f"Skipping chromosome {chr_id} (organism: {organism}) due to size < 0.05M bases")
            return None
        if chr_id.startswith('NW'):
            print(f"Skipping chromosome {chr_id} (organism: {organism}) due to NW_ prefix")
            return None

        # Memory check
        estimated_mem_gb = len(seq) / (1024 ** 3)
        if estimated_mem_gb > 31:
            print(f"Skipping chromosome {chr_id} (organism: {organism}) due to high memory requirement: {estimated_mem_gb:.2f} GB estimated")
            return None

        # Process chromosome
        exon_mask_seq, Gene_P_counter, Gene_C_counter, Exon_P_counter, Exon_C_counter, CDS_P_counter, CDS_C_counter = \
            exon_and_intron_masking(seq, df_all, chr_id)
        if CDS_P_counter == 0 or CDS_C_counter == 0:
            print(f"Skipping chromosome {chr_id} (organism: {organism}) due to zero CDS counts")
            return None

        # Initialize lists for this chromosome (kept your original data structure)
        A2T_list = []
        G2C_list = []
        P2C_number_list = [[], [], [], [], []]
        GC_countent_list_W = []
        GC_countent_list_P = [[], [], [], [], []]
        GC_countent_list_C = [[], [], [], [], []]
        len_seq_list = []
        organims_name = [organism]
        A2T_anot_list = [[], [], [], [], []]
        G2C_anot_list = [[], [], [], [], []]
        P_element_number_list = [[], [], [], [], []]
        C_element_number_list = [[], [], [], [], []]
        P_A_usage_list = [[], [], [], [], []]
        P_T_usage_list = [[], [], [], [], []]
        P_C_usage_list = [[], [], [], [], []]
        P_G_usage_list = [[], [], [], [], []]
        C_A_usage_list = [[], [], [], [], []]
        C_T_usage_list = [[], [], [], [], []]
        C_C_usage_list = [[], [], [], [], []]
        C_G_usage_list = [[], [], [], [], []]
        smaller_portion_2_Bigger = [[], [], [], [], []]
        dif2sum = [[], [], [], [], []]
        P_Portion_list = [[], [], [], [], []]
        C_portion_list = [[], [], [], [], []]
        P2C_list = [[], [], [], [], []]
        A2T_P_list = [[], [], [], [], []]  # Added
        G2C_P_list = [[], [], [], [], []]  # Added
        A2T_C_list = [[], [], [], [], []]  # Added
        G2C_C_list = [[], [], [], [], []]  # Added
        chr_id_list = [chr_id]
        All_annotation = []
        Int_annotaion = []
        UTR_annotation = []
        CDS_annotation = []
        Shared_annotation = []
        Un_annotation = []
        Sum_E = []
        Gene_P_counter_list = [Gene_P_counter]
        Gene_C_counter_list = [Gene_C_counter]
        Exon_P_counter_list = [Exon_P_counter]
        Exon_C_counter_list = [Exon_C_counter]
        CDS_P_counter_list = [CDS_P_counter]
        CDS_C_counter_list = [CDS_C_counter]

        # Nucleotide counts
        A_W = seq.count('A') + seq.count("a")
        T_W = seq.count('T') + seq.count("t")
        C_W = seq.count('C') + seq.count("c")
        G_W = seq.count("G") + seq.count("g")
        N_W = seq.count('N') + seq.count("n")
        Length = len(seq) - N_W
        denom_W = (A_W + T_W + C_W)
        if denom_W == 0:
            GC_countent_list_W.append(0)
        else:
            GC_countent_list_W.append(round((G_W + C_W) / denom_W * 100, 4))
        len_seq_list.append(round(Length / 1000000, 6))
        if T_W == 0:
            A2T_list.append(0)
        else:
            A2T_list.append(round((A_W) / (T_W), 4))
        if C_W == 0:
            G2C_list.append(0)
        else:
            G2C_list.append(round((G_W) / (C_W), 4))

        L = len(exon_mask_seq)
        df = nucleotide_analysis_for_exon(exon_mask_seq)
        # add columns "whole" and fill the same as length (this was in your original code)
        df.loc["Chr_Len"] = {"A": L, 'T': L, "C": L, "G": L}
        n_count = exon_mask_seq.count('n') + exon_mask_seq.count('N')
        df.loc["Nn_number"] = {"A": n_count, 'T': n_count, "C": n_count, "G": n_count}
        df.loc["organism"] = {"A": organism, 'T': organism, "C": organism, "G": organism}
        df.loc["chr_id"] = {"A": chr_id, 'T': chr_id, "C": chr_id, "G": chr_id}

        gene_elements = ['Int', 'UTR', 'CDS', 'Share', 'Unanot']
        counter = 0
        portions = []
        main_dict = {
            'Name': organims_name,
            'Chr_ID': chr_id_list,
            'Length(M)': len_seq_list,
            'GC_W': GC_countent_list_W,
            'A2T_W': A2T_list,
            'G2C_W': G2C_list,
        }

        for e_number, element in enumerate(gene_elements):
            if element == "Share" or element == "Unanot":
                Ap, Tp, Cp, Gp = df.iloc[counter + 0]
                Ap, Tp, Cp, Gp = Ap / 2, Tp / 2, Cp / 2, Gp / 2
                Ac, Tc, Cc, Gc = df.iloc[counter + 0]
                Ac, Tc, Cc, Gc = Ac / 2, Tc / 2, Cc / 2, Gc / 2
                counter += 1
            else:
                Ap, Tp, Cp, Gp = df.iloc[counter + 0]
                Ac, Tc, Cc, Gc = df.iloc[counter + 1]
                counter += 2
            try:
                # you had two lines setting u — left them in place but they don't affect outputs
                u = 1 / Ap + 1 / Tp + 1 / Cp + 1 / Gp
                u = 1 / Ac + 1 / Tc + 1 / Cc + 1 / Gc
                portions.append(sum([Ap, Tp, Cp, Gp, Ac, Tc, Cc, Gc]) / Length * 100)

                P_portion = round(Ap + Tp + Cp + Gp, 4)
                if P_portion == 0:
                    P_A_usage_list[e_number].append(0)
                    P_T_usage_list[e_number].append(0)
                    P_C_usage_list[e_number].append(0)
                    P_G_usage_list[e_number].append(0)
                    P_Portion_list[e_number].append(0)
                    A2T_P_list[e_number].append(0)
                    G2C_P_list[e_number].append(0)
                    GC_countent_list_P[e_number].append(0)
                else:
                    P_A_usage_list[e_number].append(round(Ap / P_portion * 100, 4))
                    P_T_usage_list[e_number].append(round(Tp / P_portion * 100, 4))
                    P_C_usage_list[e_number].append(round(Cp / P_portion * 100, 4))
                    P_G_usage_list[e_number].append(round(Gp / P_portion * 100, 4))
                    P_Portion_list[e_number].append(round(P_portion / Length * 100, 4))
                    if Tp == 0:
                        A2T_P_list[e_number].append(0)
                    else:
                        A2T_P_list[e_number].append(round(Ap / Tp, 4))
                    if Cp == 0:
                        G2C_P_list[e_number].append(0)
                    else:
                        G2C_P_list[e_number].append(round(Gp / Cp, 4))
                    P_element_number_list[e_number].append(CDS_P_counter)
                    GC_countent_list_P[e_number].append(round((Cp + Gp) / P_portion * 100, 4))

                C_portion = round(Ac + Tc + Cc + Gc, 4)
                if C_portion == 0:
                    C_A_usage_list[e_number].append(0)
                    C_T_usage_list[e_number].append(0)
                    C_C_usage_list[e_number].append(0)
                    C_G_usage_list[e_number].append(0)
                    C_portion_list[e_number].append(0)
                    A2T_C_list[e_number].append(0)
                    G2C_C_list[e_number].append(0)
                    GC_countent_list_C[e_number].append(0)
                else:
                    C_A_usage_list[e_number].append(round(Ac / C_portion * 100, 4))
                    C_T_usage_list[e_number].append(round(Tc / C_portion * 100, 4))
                    C_C_usage_list[e_number].append(round(Cc / C_portion * 100, 4))
                    C_G_usage_list[e_number].append(round(Gc / C_portion * 100, 4))
                    C_portion_list[e_number].append(round(C_portion / Length * 100, 4))
                    if Tc == 0:
                        A2T_C_list[e_number].append(0)
                    else:
                        A2T_C_list[e_number].append(round(Ac / Tc, 4))
                    if Cc == 0:
                        G2C_C_list[e_number].append(0)
                    else:
                        G2C_C_list[e_number].append(round(Gc / Cc, 4))
                    C_element_number_list[e_number].append(CDS_C_counter)
                    GC_countent_list_C[e_number].append(round((Cc + Gc) / C_portion * 100, 4))

                A2T_anot_list[e_number].append(round((Ap + Ac) / (Tp + Tc) if (Tp + Tc) != 0 else 0, 4))
                G2C_anot_list[e_number].append(round((Gp + Gc) / (Cp + Cc) if (Cp + Cc) != 0 else 0, 4))

                P2C_list[e_number].append(round(P_portion / C_portion, 4) if C_portion != 0 else 0)
                P2C_number_list[e_number].append(round(CDS_P_counter / CDS_C_counter, 4) if CDS_C_counter != 0 else 0)
                PandC = [P_portion, C_portion]
                PandC.sort()
                if PandC[1] == 0:
                    smaller_portion_2_Bigger[e_number].append(0)
                else:
                    smaller_portion_2_Bigger[e_number].append(round(PandC[0] / PandC[1], 4))
                dif2sum[e_number].append(round(abs(PandC[0] - PandC[1]) / ((PandC[0] + PandC[1]) * 2) * 100, 4) if (PandC[0] + PandC[1]) != 0 else 0)

            except ZeroDivisionError:
                portions.append(0)
                P_portion = round(Ap + Tp + Cp + Gp, 4)
                P_A_usage_list[e_number].append(0)
                P_T_usage_list[e_number].append(0)
                P_C_usage_list[e_number].append(0)
                P_G_usage_list[e_number].append(0)
                P_Portion_list[e_number].append(0)
                A2T_P_list[e_number].append(0)
                G2C_P_list[e_number].append(0)
                P_element_number_list[e_number].append(0)
                GC_countent_list_P[e_number].append(0)

                C_portion = round(Ac + Tc + Cc + Gc, 4)
                C_A_usage_list[e_number].append(0)
                C_T_usage_list[e_number].append(0)
                C_C_usage_list[e_number].append(0)
                C_G_usage_list[e_number].append(0)
                C_portion_list[e_number].append(0)
                A2T_C_list[e_number].append(0)
                G2C_C_list[e_number].append(0)
                C_element_number_list[e_number].append(0)
                GC_countent_list_C[e_number].append(0)

                A2T_anot_list[e_number].append(0)
                G2C_anot_list[e_number].append(0)

                P2C_list[e_number].append(0)
                P2C_number_list[e_number].append(0)
                PandC = [P_portion, C_portion]
                PandC.sort()
                smaller_portion_2_Bigger[e_number].append(0)
                dif2sum[e_number].append(0)

            vaiavle_dic = {
                f'{element}_GC_P': GC_countent_list_P[e_number],
                f'{element}_GC_C': GC_countent_list_C[e_number],
                f'{element}_P%': P_Portion_list[e_number],
                f'{element}_C%': C_portion_list[e_number],
                f'{element}_P2C': P2C_list[e_number],
                f'{element}_S2B': smaller_portion_2_Bigger[e_number],
                f'{element}_dif2sum': dif2sum[e_number],
                f'{element}_A2T_ano': A2T_anot_list[e_number],
                f'{element}_G2C_ano': G2C_anot_list[e_number],
                f'{element}_P_E_num': P_element_number_list[e_number],
                f'{element}_C_E_num': C_element_number_list[e_number],
                f'{element}_P_A_ano': P_A_usage_list[e_number],
                f'{element}_P_T_ano': P_T_usage_list[e_number],
                f'{element}_P_C_ano': P_C_usage_list[e_number],
                f'{element}_P_G_ano': P_G_usage_list[e_number],
                f'{element}_P_A2T': A2T_P_list[e_number],
                f'{element}_P_G2C': G2C_P_list[e_number],
                f'{element}_C_A_ano': C_A_usage_list[e_number],
                f'{element}_C_T_ano': C_T_usage_list[e_number],
                f'{element}_C_C_ano': C_C_usage_list[e_number],
                f'{element}_C_G_ano': C_G_usage_list[e_number],
                f'{element}_A2T_C': A2T_C_list[e_number],
                f'{element}_G2C_C': G2C_C_list[e_number]
            }

            main_dict = main_dict | vaiavle_dic

        All_annotation.append(round(sum(portions[:4]), 4))
        Int_annotaion.append(round(portions[0], 4))
        UTR_annotation.append(round(portions[1], 4))
        CDS_annotation.append(round(portions[2], 4))
        Shared_annotation.append(round(portions[3], 4))
        Un_annotation.append(round(portions[4], 4))
        Sum_E.append(round(sum(portions), 4))

        try:
            my_table = pd.DataFrame(main_dict)
        except ValueError as e:
            print(f"Error creating DataFrame for {organism} {chr_id}: {e}")
            for i in main_dict:
                print('hi', i, len(main_dict[i]))
            return None

        my_table.insert(4, 'Wh_ano', All_annotation)
        my_table.insert(5, 'Int_ano', Int_annotaion)
        my_table.insert(6, 'UTR_ano', UTR_annotation)
        my_table.insert(7, 'CDS_ano', CDS_annotation)
        my_table.insert(8, 'share_ano', Shared_annotation)
        my_table.insert(9, 'Un_ano', Un_annotation)
        my_table.insert(10, 'sum', Sum_E)
        my_table['Gene_P_count'] = Gene_P_counter_list
        my_table['Gene_C_count'] = Gene_C_counter_list
        my_table['Exon_P_count'] = Exon_P_counter_list
        my_table['Exon_C_count'] = Exon_C_counter_list
        my_table['CDS_P_count'] = CDS_P_counter_list
        my_table['CDS_C_count'] = CDS_C_counter_list

        # Print header and row as in original
        column_width = 15
        header = " ".join(col.ljust(column_width) for col in my_table.columns)
        last_row_values = " ".join(str(val).ljust(column_width) for val in my_table.iloc[-1].values)
        print(header)
        print("-" * len(header))
        print(final_number, chr_counter, last_row_values)

        return my_table

    except Exception as e:
        print(f"Error in chromosome {chr_id} (organism: {organism}): {e}")
        return None

# ----------------- Helper Function -----------------
def load_chromosome_batch(all_genomes_folder_path: str, remaining_organisms: List[str], batch_size: int, gff_cache: dict) -> List[Tuple[int, str, str]]:
    """
    Load up to batch_size chromosomes from remaining organisms.
    Returns: List of (size, organism, chr_id) tuples for valid chromosomes.
    Caches GFF DataFrames in gff_cache to avoid reloading.
    """
    chromosomes = []
    current_org_idx = 0

    while len(chromosomes) < batch_size and current_org_idx < len(remaining_organisms):
        org = remaining_organisms[current_org_idx]
        try:
            files = os.listdir(os.path.join(all_genomes_folder_path, org))
            fna = [i for i in files if i.endswith('.fna.gz') and not i.startswith('cds')]
            if not fna:
                current_org_idx += 1
                continue
            fasta_path = os.path.join(all_genomes_folder_path, org, fna[0])

            # Load GFF if not cached
            if org not in gff_cache:
                gff = [k for k in files if k.endswith('.gff.gz')]
                if not gff:
                    current_org_idx += 1
                    continue
                gff_path = os.path.join(all_genomes_folder_path, org, gff[0])
                if gff_path.endswith('.gz'):
                    with gzip.open(gff_path, 'rt') as gff_file:
                        gff_cache[org] = pd.read_csv(gff_file, sep='\t', comment='#')
                else:
                    with open(gff_path, 'r') as gff_file:
                        gff_cache[org] = pd.read_csv(gff_file, sep='\t', comment='#')

            # Scan chromosomes
            with gzip.open(fasta_path, 'rt') as fasta_file:
                for seq_record in SeqIO.parse(fasta_file, "fasta"):
                    if 'mitochondrion' in seq_record.description:
                        continue
                    if "chromosome" not in seq_record.description:
                        continue
                    if seq_record.id.startswith('NW'):
                        continue
                    seq_len = len(seq_record.seq)
                    if seq_len > 300000000 or seq_len / 1000000 < 0.05:
                        continue
                    chromosomes.append((seq_len, org, seq_record.id))
                    if len(chromosomes) >= batch_size:
                        break
            current_org_idx += 1
        except Exception as e:
            print(f"Error scanning {org}: {e}")
            current_org_idx += 1
            continue

    return chromosomes

# ----------------- Parallel Controller -----------------
def gene_element_annotation_in_plus_and_minis(all_genomes_folder_path: str, organ_group: str):
    print('Analysis started')
    start_time = time.time()
    current_time = datetime.now()
    start_run = current_time.strftime("%Y-%m-%d_%H-%M-%S")
    print(start_run)

    lablename = 'refrence'
    ram_limit = RAM_LIMIT_GB * (1024 ** 3) * USABLE_FRACTION  # ~22.3GB

    # Get list of organisms
    all_organisms = [org for org in os.listdir(all_genomes_folder_path)
                     if not (org.endswith(".json") or org.endswith(".jsonl"))]
    organism_number = len(all_organisms)
    print(f"Found {organism_number} organisms to process (estimating {BATCH_SIZE} chromosomes per batch)")

    # Queues and trackers
    ready_queue = []  # Priority queue: [(-size, org, chr_id)]
    remaining_organisms = all_organisms[:]  # Copy of unprocessed organisms
    gff_cache = {}  # Cache GFF DataFrames to avoid reloading
    results = []
    active_tasks = []  # List of (future, est_mem, org, chr_id)
    chr_counter = 0  # Global chromosome counter

    def load_next_batch():
        """Load up to BATCH_SIZE chromosomes into ready_queue."""
        chromosomes = load_chromosome_batch(all_genomes_folder_path, remaining_organisms, BATCH_SIZE, gff_cache)
        if chromosomes:  # Check if batch is non-empty
            for size, org, chr_id in sorted(chromosomes, key=lambda x: x[0], reverse=True):
                heapq.heappush(ready_queue, (-size, org, chr_id))  # Negative for max-heap
            print(f"Batch loaded: {len(chromosomes)} chromosomes, largest: {chromosomes[0][0]:,} bases")
        else:
            print("No chromosomes loaded in this batch")

    with ProcessPoolExecutor(max_workers=MAX_CORES) as executor:
        # Load initial batch
        load_next_batch()

        while ready_queue or active_tasks:
            # Submit next task if possible and ready_queue has items
            if ready_queue:
                neg_size, org, chr_id = heapq.heappop(ready_queue)
                size = -neg_size
                est_mem = (size * BYTES_PER_BASE) + PROCESS_OVERHEAD
                current_free_gb = psutil.virtual_memory().available / (1024 ** 3)
                current_estimated_used = sum(mem for _, mem, _, _ in active_tasks)

                # Check if we can submit
                if (current_free_gb * BUFFER_FACTOR > est_mem / (1024 ** 3) and
                        current_estimated_used + est_mem < ram_limit):
                    # Load sequence for this chromosome
                    fasta_path = os.path.join(all_genomes_folder_path, org, [f for f in os.listdir(os.path.join(all_genomes_folder_path, org)) if f.endswith('.fna.gz') and not f.startswith('cds')][0])
                    with gzip.open(fasta_path, 'rt') as fasta_file:
                        seq = None
                        for seq_record in SeqIO.parse(fasta_file, "fasta"):
                            if seq_record.id == chr_id:
                                seq = str(seq_record.seq)
                                break
                        if seq is None:
                            print(f"Chromosome {chr_id} not found in {org}")
                            continue

                    future = executor.submit(
                        process_one_chromosome, organism_number, chr_counter, org, chr_id, seq, gff_cache[org], all_genomes_folder_path
                    )
                    active_tasks.append((future, est_mem, org, chr_id))
                    #print(f"Submitted {org} {chr_id} (size: {size:,}, est. mem: {est_mem / (1024**2):.1f} MB). Active: {len(active_tasks)}, Queue: {len(ready_queue)}")
                    chr_counter += 1
                else:
                    # Can't submit yet; push back to queue
                    heapq.heappush(ready_queue, (neg_size, org, chr_id))
                    #print(f"Waiting for RAM... Free: {current_free_gb:.1f}GB, Est. used: {current_estimated_used / (1024**3):.1f}GB")

            # Load next batch if queue is low
            if len(ready_queue) < BATCH_SIZE // 2 and remaining_organisms:
                load_next_batch()

            # Check for completions
            if active_tasks:
                # WAIT for at least one future to complete (non-blocking with timeout)
                done, _ = wait(
                    [f for f, _, _, _ in active_tasks], timeout=0.5, return_when=FIRST_COMPLETED
                )
                for future in done:
                    for i, (f, mem, org, chr_id) in enumerate(active_tasks):
                        if f == future:
                            try:
                                res = future.result()
                                if isinstance(res, pd.DataFrame):
                                    results.append(res)
                                print(f"Completed {org} {chr_id} (freed {mem / (1024**2):.1f} MB). Active: {len(active_tasks)-1}")
                            except Exception as exc:
                                print(f"{org} {chr_id} failed: {exc}")
                            del active_tasks[i]
                            break
            else:
                time.sleep(0.1)  # Idle wait

    # Clean up GFF cache
    gff_cache.clear()

    end_time = time.time()
    print('Process time:', round((end_time - start_time) / 60, 2), 'minutes')

    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    results_folder = f"{organ_group}_{lablename}_Int_UTR_CDS_{timestamp}"
    os.makedirs(results_folder, exist_ok=True)

    if results:
        try:
            final_df = pd.concat(results, ignore_index=True)
            output_filename = os.path.join(results_folder, f"{organ_group}_results_{timestamp}.xlsx")
            final_df.to_excel(output_filename, index=False)
            print(f"Results saved to: {output_filename}")
        except Exception as e:
            print("Error saving final Excel:", e)
    else:
        print("No results to save.")

    print('Start running:', start_run)
    print('End running:', timestamp)
    print('Process finished')

# ----------------- Run -----------------
if __name__ == "__main__":
    one_Trabite = r"/media/su/Bioinf/Genomes_2025_3/Genomes_2025_3_7/chordata_annotated_complete_chr_genome_gff_cds_2025/ncbi_dataset/data"
    gene_element_annotation_in_plus_and_minis(one_Trabite, 'Chordata')