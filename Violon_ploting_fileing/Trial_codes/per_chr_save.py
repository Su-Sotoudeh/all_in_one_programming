#!/usr/bin/env python3
"""
Parallel genome analysis — row-by-row streaming to single CSV with resume.
Produces 134 columns per chromosome (exact order requested).
"""

import os
import gzip
import time
import math
import multiprocessing
from multiprocessing import Manager, Process
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
import pandas as pd
from Bio import SeqIO

# ---------------- USER CONFIG ----------------
ONE_TRABITE = r"/media/su/Bioinf/Genomes_2025_3/Genomes_2025_3_7/Bacteria_2025_2_1/bacteria_reference_complet_genomes/ncbi_dataset/data"
NUM_WORKERS = None   # None -> uses all CPUs
RESULTS_FILE = f"genome_results_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.csv"

# ---------------- TRANSLATION TABLES (bytearray.translate) ----------------
def build_translation_tables():
    t = {}
    # gene plus
    gene_plus = bytearray(range(256))
    gene_plus[ord('a')] = ord('Q')
    gene_plus[ord('t')] = ord('W')
    gene_plus[ord('c')] = ord('E')
    gene_plus[ord('g')] = ord('R')
    gene_plus[ord('q')] = ord('A')
    gene_plus[ord('w')] = ord('T')
    gene_plus[ord('e')] = ord('C')
    gene_plus[ord('r')] = ord('G')
    # gene minus
    gene_minus = bytearray(range(256))
    gene_minus[ord('a')] = ord('q')
    gene_minus[ord('t')] = ord('w')
    gene_minus[ord('c')] = ord('e')
    gene_minus[ord('g')] = ord('r')
    gene_minus[ord('Q')] = ord('A')
    gene_minus[ord('W')] = ord('T')
    gene_minus[ord('E')] = ord('C')
    gene_minus[ord('R')] = ord('G')

    # exon plus
    exon_plus = bytearray(range(256))
    exon_plus[ord('Q')] = ord('S')
    exon_plus[ord('W')] = ord('D')
    exon_plus[ord('E')] = ord('F')
    exon_plus[ord('R')] = ord('H')
    exon_plus[ord('q')] = ord('A')
    exon_plus[ord('w')] = ord('T')
    exon_plus[ord('e')] = ord('C')
    exon_plus[ord('r')] = ord('G')
    # exon minus
    exon_minus = bytearray(range(256))
    exon_minus[ord('q')] = ord('s')
    exon_minus[ord('w')] = ord('d')
    exon_minus[ord('e')] = ord('f')
    exon_minus[ord('r')] = ord('h')
    exon_minus[ord('Q')] = ord('A')
    exon_minus[ord('W')] = ord('T')
    exon_minus[ord('E')] = ord('C')
    exon_minus[ord('R')] = ord('G')

    # cds plus
    cds_plus = bytearray(range(256))
    cds_plus[ord('S')] = ord('Z')
    cds_plus[ord('D')] = ord('X')
    cds_plus[ord('F')] = ord('V')
    cds_plus[ord('H')] = ord('B')
    cds_plus[ord('z')] = ord('A')
    cds_plus[ord('x')] = ord('T')
    cds_plus[ord('v')] = ord('C')
    cds_plus[ord('b')] = ord('G')
    # cds minus
    cds_minus = bytearray(range(256))
    cds_minus[ord('s')] = ord('z')
    cds_minus[ord('d')] = ord('x')
    cds_minus[ord('f')] = ord('v')
    cds_minus[ord('h')] = ord('b')
    cds_minus[ord('Z')] = ord('A')
    cds_minus[ord('X')] = ord('T')
    cds_minus[ord('V')] = ord('C')
    cds_minus[ord('B')] = ord('G')

    t["gene_plus"] = gene_plus
    t["gene_minus"] = gene_minus
    t["exon_plus"] = exon_plus
    t["exon_minus"] = exon_minus
    t["cds_plus"] = cds_plus
    t["cds_minus"] = cds_minus
    return t

TABLES = build_translation_tables()

# ----------------- column order (exact) -----------------
COLUMN_ORDER = [
    "Name","Chr_ID","Length(M)","GC_W","Wh_ano","Int_ano","UTR_ano","CDS_ano","share_ano","Un_ano","sum",
    "A2T_W","G2C_W",
    # Int fields
    "Int_GC_P","Int_GC_C","Int_P%","Int_C%","Int_P2C","Int_S2B","Int_dif2sum","Int_A2T_ano","Int_G2C_ano",
    "Int_P_E_num","Int_C_E_num",
    "Int_P_A_ano","Int_P_T_ano","Int_P_C_ano","Int_P_G_ano","Int_P_A2T","Int_P_G2C",
    "Int_C_A_ano","Int_C_T_ano","Int_C_C_ano","Int_C_G_ano","Int_A2T_C","Int_G2C_C",
    # UTR
    "UTR_GC_P","UTR_GC_C","UTR_P%","UTR_C%","UTR_P2C","UTR_S2B","UTR_dif2sum","UTR_A2T_ano","UTR_G2C_ano",
    "UTR_P_E_num","UTR_C_E_num",
    "UTR_P_A_ano","UTR_P_T_ano","UTR_P_C_ano","UTR_P_G_ano","UTR_P_A2T","UTR_P_G2C",
    "UTR_C_A_ano","UTR_C_T_ano","UTR_C_C_ano","UTR_C_G_ano","UTR_A2T_C","UTR_G2C_C",
    # CDS
    "CDS_GC_P","CDS_GC_C","CDS_P%","CDS_C%","CDS_P2C","CDS_S2B","CDS_dif2sum","CDS_A2T_ano","CDS_G2C_ano",
    "CDS_P_E_num","CDS_C_E_num",
    "CDS_P_A_ano","CDS_P_T_ano","CDS_P_C_ano","CDS_P_G_ano","CDS_P_A2T","CDS_P_G2C",
    "CDS_C_A_ano","CDS_C_T_ano","CDS_C_C_ano","CDS_C_G_ano","CDS_A2T_C","CDS_G2C_C",
    # Share
    "Share_GC_P","Share_GC_C","Share_P%","Share_C%","Share_P2C","Share_S2B","Share_dif2sum","Share_A2T_ano","Share_G2C_ano",
    "Share_P_E_num","Share_C_E_num",
    "Share_P_A_ano","Share_P_T_ano","Share_P_C_ano","Share_P_G_ano","Share_P_A2T","Share_P_G2C",
    "Share_C_A_ano","Share_C_T_ano","Share_C_C_ano","Share_C_G_ano","Share_A2T_C","Share_G2C_C",
    # Unannotated
    "Unanot_GC_P","Unanot_GC_C","Unanot_P%","Unanot_C%","Unanot_P2C","Unanot_S2B","Unanot_dif2sum","Unanot_A2T_ano","Unanot_G2C_ano",
    "Unanot_P_E_num","Unanot_C_E_num",
    "Unanot_P_A_ano","Unanot_P_T_ano","Unanot_P_C_ano","Unanot_P_G_ano","Unanot_P_A2T","Unanot_P_G2C",
    "Unanot_C_A_ano","Unanot_C_T_ano","Unanot_C_C_ano","Unanot_C_G_ano","Unanot_A2T_C","Unanot_G2C_C",
    # counts
    "Gene_P_count","Gene_C_count","Exon_P_count","Exon_C_count","CDS_P_count","CDS_C_count"
]

# ----------------- helper functions -----------------
def exon_and_intron_masking(sequence, gff_path, chr_id):
    """Fast in-place annotate using bytearray.translate; returns annotated string and counters."""
    seq = bytearray(sequence.lower(), 'ascii')  # 1 byte per base
    Intron_P = Intron_C = Exon_P = Exon_C = CDS_P = CDS_C = 0

    # read gff
    if gff_path.endswith('.gz'):
        gff_f = gzip.open(gff_path, 'rt')
    else:
        gff_f = open(gff_path, 'rt')
    try:
        df = pd.read_csv(gff_f, sep='\t', comment='#', header=None,
                         names=['seqid','source','type','start','end','score','strand','phase','attributes'])
    finally:
        gff_f.close()

    ch_id_split = chr_id.split('.')[0]

    # gene entries
    df_genes = df[(df['seqid'].astype(str).str.startswith(ch_id_split)) & (df['type'].astype(str).str.startswith('gene'))]
    for _, row in df_genes.iterrows():
        start = int(row['start']) - 1
        end = int(row['end'])
        if end <= start:
            continue
        if str(row['strand']) == '+':
            seq[start:end] = seq[start:end].translate(TABLES["gene_plus"])
            Intron_P += 1
        else:
            seq[start:end] = seq[start:end].translate(TABLES["gene_minus"])
            Intron_C += 1

    # exons
    df_exons = df[(df['seqid'].astype(str).str.startswith(ch_id_split)) & (df['type'].astype(str).str.startswith('exon'))]
    for _, row in df_exons.iterrows():
        start = int(row['start']) - 1
        end = int(row['end'])
        if end <= start:
            continue
        if str(row['strand']) == '+':
            seq[start:end] = seq[start:end].translate(TABLES["exon_plus"])
            Exon_P += 1
        else:
            seq[start:end] = seq[start:end].translate(TABLES["exon_minus"])
            Exon_C += 1

    # CDS
    df_cds = df[(df['seqid'].astype(str).str.startswith(ch_id_split)) & (df['type'].astype(str).str.startswith('CDS'))]
    for _, row in df_cds.iterrows():
        start = int(row['start']) - 1
        end = int(row['end'])
        if end <= start:
            continue
        if str(row['strand']) == '+':
            seq[start:end] = seq[start:end].translate(TABLES["cds_plus"])
            CDS_P += 1
        else:
            seq[start:end] = seq[start:end].translate(TABLES["cds_minus"])
            CDS_C += 1

    return seq.decode('ascii'), Intron_P, Intron_C, Exon_P, Exon_C, CDS_P, CDS_C

def nucleotide_analysis_for_exon(sequence_plus):
    """Return df with rows: +entron,-entron,+UTR,-UTR,+CDS,-CDS,Shared,Unannotated"""
    A = [sequence_plus.count(ch) for ch in ['Q','q','S','s','Z','z','A','a']]
    T = [sequence_plus.count(ch) for ch in ['W','w','D','d','X','x','T','t']]
    C = [sequence_plus.count(ch) for ch in ['E','e','F','f','V','v','C','c']]
    G = [sequence_plus.count(ch) for ch in ['R','r','H','h','B','b','G','g']]
    df = pd.DataFrame({'A':A,'T':T,'C':C,'G':G}, index=['+entron','-entron','+UTR','-UTR','+CDS','-CDS','Shared','Unannotated'])
    return df

# writer process — single writer to avoid file corruption
def writer_proc(queue, results_file, header_columns, existing_chr):
    """
    Reads row-dicts from queue and appends them to results_file in CSV.
    Skips writing if Chr_ID already in existing_chr or already_written set.
    """
    already_written = set(existing_chr) if existing_chr else set()
    first_write = not os.path.exists(results_file)
    # Open file handle in append mode, but we'll rely on pandas to_csv for header behaviour
    while True:
        item = queue.get()
        if item is None:
            # sentinel -> finish
            break
        # item expected to be dict
        chr_id = str(item.get('Chr_ID',''))
        if chr_id in already_written:
            # skip duplicate row
            continue
        # Ensure dict has all columns in order; missing -> 0
        row = {col: item.get(col, 0) for col in header_columns}
        df_row = pd.DataFrame([row])
        try:
            df_row.to_csv(results_file, mode='a', header=first_write, index=False)
            first_write = False
            already_written.add(chr_id)
        except Exception as e:
            print("Writer error saving row:", e)
    # writer finishing
    return

# Build a single-row dictionary of 134 columns for a chromosome
def build_row_from_counts(organism, ch_id, seq, exon_mask_seq, df_ann, gene_counts):
    """
    df_ann: DataFrame from nucleotide_analysis_for_exon
    gene_counts: tuple (Gene_P, Gene_C, Exon_P, Exon_C, CDS_P, CDS_C)
    """
    Gene_P, Gene_C, Exon_P, Exon_C, CDS_P, CDS_C = gene_counts
    # base counts on whole seq
    A_W = seq.count('A') + seq.count('a')
    T_W = seq.count('T') + seq.count('t')
    C_W = seq.count('C') + seq.count('c')
    G_W = seq.count('G') + seq.count('g')
    N_W = seq.count('N') + seq.count('n')
    Length = len(seq) - N_W if (len(seq) - N_W) > 0 else 1  # avoid zero

    # compute portions per element
    gene_elements = ['Int','UTR','CDS','Share','Unanot']
    counter = 0
    portions = []
    element_values = {}  # store per-element metrics

    for e_number, element in enumerate(gene_elements):
        if element == "Share" or element == "Unanot":
            # share/unannotated use single row halved for P and C
            Ap, Tp, Cp, Gp = df_ann.iloc[counter + 0].tolist()
            Ap, Tp, Cp, Gp = Ap / 2.0, Tp / 2.0, Cp / 2.0, Gp / 2.0
            Ac, Tc, Cc, Gc = df_ann.iloc[counter + 0].tolist()
            Ac, Tc, Cc, Gc = Ac / 2.0, Tc / 2.0, Cc / 2.0, Gc / 2.0
            counter += 1
        else:
            Ap, Tp, Cp, Gp = df_ann.iloc[counter + 0].tolist()
            Ac, Tc, Cc, Gc = df_ann.iloc[counter + 1].tolist()
            counter += 2

        # ensure floats
        Ap, Tp, Cp, Gp = float(Ap), float(Tp), float(Cp), float(Gp)
        Ac, Tc, Cc, Gc = float(Ac), float(Tc), float(Cc), float(Gc)

        P_portion = Ap + Tp + Cp + Gp
        C_portion = Ac + Tc + Cc + Gc

        # portion percent of chromosome occupied by this annotation (P+C)/Length*100
        try:
            portion_percent = (P_portion + C_portion) / Length * 100.0
        except Exception:
            portion_percent = 0.0
        portions.append(portion_percent)

        # safe computations with zero checks
        # P usages %
        P_A_ano = round(Ap / P_portion * 100, 4) if P_portion != 0 else 0
        P_T_ano = round(Tp / P_portion * 100, 4) if P_portion != 0 else 0
        P_C_ano = round(Cp / P_portion * 100, 4) if P_portion != 0 else 0
        P_G_ano = round(Gp / P_portion * 100, 4) if P_portion != 0 else 0
        P_percent = round(P_portion / Length * 100, 4) if Length != 0 else 0
        # C usages %
        C_A_ano = round(Ac / C_portion * 100, 4) if C_portion != 0 else 0
        C_T_ano = round(Tc / C_portion * 100, 4) if C_portion != 0 else 0
        C_C_ano = round(Cc / C_portion * 100, 4) if C_portion != 0 else 0
        C_G_ano = round(Gc / C_portion * 100, 4) if C_portion != 0 else 0
        C_percent = round(C_portion / Length * 100, 4) if Length != 0 else 0

        # GC content inside P and C
        GC_P = round(((Cp + Gp) / P_portion) * 100, 4) if P_portion != 0 else 0
        GC_C = round(((Cc + Gc) / C_portion) * 100, 4) if C_portion != 0 else 0

        # A2T, G2C ratios inside P and C
        P_A2T = round(Ap / Tp, 4) if Tp != 0 else 0
        P_G2C = round(Gp / Cp, 4) if Cp != 0 else 0
        C_A2T = round(Ac / Tc, 4) if Tc != 0 else 0
        C_G2C = round(Gc / Cc, 4) if Cc != 0 else 0

        # A2T_ano and G2C_ano across P+C
        A2T_ano = round((Ap + Ac) / (Tp + Tc), 4) if (Tp + Tc) != 0 else 0
        G2C_ano = round((Gp + Gc) / (Cp + Cc), 4) if (Cp + Cc) != 0 else 0

        # P2C, counts, S2B, dif2sum
        P2C = round(P_portion / C_portion, 4) if C_portion != 0 else 0
        P_E_num = CDS_P  # original code used CDS_P_counter for these counts
        C_E_num = CDS_C
        if max(P_portion, C_portion) != 0:
            S2B = round(min(P_portion, C_portion) / max(P_portion, C_portion), 4)
        else:
            S2B = 0
        if (P_portion + C_portion) != 0:
            dif2sum = round(abs(P_portion - C_portion) / ((P_portion + C_portion) * 2) * 100, 4)
        else:
            dif2sum = 0

        # store values in element_values with appropriate prefixes
        prefix = element
        element_values[f"{prefix}_GC_P"] = GC_P
        element_values[f"{prefix}_GC_C"] = GC_C
        element_values[f"{prefix}_P%"] = P_percent
        element_values[f"{prefix}_C%"] = C_percent
        element_values[f"{prefix}_P2C"] = P2C
        element_values[f"{prefix}_S2B"] = S2B
        element_values[f"{prefix}_dif2sum"] = dif2sum
        element_values[f"{prefix}_A2T_ano"] = A2T_ano
        element_values[f"{prefix}_G2C_ano"] = G2C_ano
        element_values[f"{prefix}_P_E_num"] = P_E_num
        element_values[f"{prefix}_C_E_num"] = C_E_num
        element_values[f"{prefix}_P_A_ano"] = P_A_ano
        element_values[f"{prefix}_P_T_ano"] = P_T_ano
        element_values[f"{prefix}_P_C_ano"] = P_C_ano
        element_values[f"{prefix}_P_G_ano"] = P_G_ano
        element_values[f"{prefix}_P_A2T"] = P_A2T
        element_values[f"{prefix}_P_G2C"] = P_G2C
        element_values[f"{prefix}_C_A_ano"] = C_A_ano
        element_values[f"{prefix}_C_T_ano"] = C_T_ano
        element_values[f"{prefix}_C_C_ano"] = C_C_ano
        element_values[f"{prefix}_C_G_ano"] = C_G_ano
        element_values[f"{prefix}_A2T_C"] = C_A2T
        element_values[f"{prefix}_G2C_C"] = C_G2C

    # summarize top-level fields
    Wh_ano = round(sum(portions[:4]), 4) if len(portions) >= 4 else round(sum(portions), 4)
    Int_ano = round(portions[0], 4) if len(portions) > 0 else 0
    UTR_ano = round(portions[1], 4) if len(portions) > 1 else 0
    CDS_ano = round(portions[2], 4) if len(portions) > 2 else 0
    Share_ano = round(portions[3], 4) if len(portions) > 3 else 0
    Un_ano = round(portions[4], 4) if len(portions) > 4 else 0
    Sum_E = round(sum(portions), 4) if portions else 0

    # whole-genome ratios
    A2T_W = round((A_W) / (T_W), 4) if T_W != 0 else 0
    G2C_W = round((G_W) / (C_W), 4) if C_W != 0 else 0
    GC_W = round((G_W + C_W) / (A_W + T_W + C_W) * 100, 4) if (A_W + T_W + C_W) != 0 else 0

    # build final row dict following COLUMN_ORDER
    row = {}
    row["Name"] = organism
    row["Chr_ID"] = ch_id
    row["Length(M)"] = round(Length / 1e6, 6)
    row["GC_W"] = GC_W
    row["Wh_ano"] = Wh_ano
    row["Int_ano"] = Int_ano
    row["UTR_ano"] = UTR_ano
    row["CDS_ano"] = CDS_ano
    row["share_ano"] = Share_ano
    row["Un_ano"] = Un_ano
    row["sum"] = Sum_E
    row["A2T_W"] = A2T_W
    row["G2C_W"] = G2C_W

    # fill per-element fields in the exact column names
    for prefix in ["Int","UTR","CDS","Share","Unanot"]:
        # note: in your example Unannotated used prefix "Unanot"
        # element_values used key like "Unanot_GC_P" for Unanot -> ensure mapping
        mapped_prefix = prefix if prefix != "Unanot" else "Unanot"
        # some keys in element_values may use "Unanot_*"
        for key in [f"{prefix}_GC_P", f"{prefix}_GC_C", f"{prefix}_P%", f"{prefix}_C%", f[f"{prefix}_P2C" if False else f"{prefix}_P2C"] if False else f"{prefix}_P2C"]:
            pass  # placeholder to keep consistent style (we will directly use COLUMN_ORDER below)

    # finalize mapping by iterating columns
    for col in COLUMN_ORDER:
        if col in row:
            continue
        if col in element_values:
            row[col] = element_values[col]
        else:
            # if not present (should not happen), set 0
            row[col] = 0

    # counts
    row["Gene_P_count"] = Gene_P
    row["Gene_C_count"] = Gene_C
    row["Exon_P_count"] = Exon_P
    row["Exon_C_count"] = Exon_C
    row["CDS_P_count"] = CDS_P
    row["CDS_C_count"] = CDS_C

    # Round all floats to 4 decimals except some that represent counts
    for k, v in row.items():
        # keep integer counts as int
        if k.endswith("_count") or k in ("Gene_P_count","Gene_C_count","Exon_P_count","Exon_C_count","CDS_P_count","CDS_C_count"):
            continue
        try:
            if isinstance(v, float):
                row[k] = round(v, 4)
        except Exception:
            pass

    return row

# ----------------- worker per organism -----------------
def process_one_organism(organism, all_genomes_folder_path, queue, existing_completed):
    """
    Process one organism; for each chromosome produce a row dict and put into queue.
    existing_completed: set of Chr_ID (strings) already in CSV at start time (resume).
    """
    try:
        org_dir = os.path.join(all_genomes_folder_path, organism)
        files = os.listdir(org_dir)
        fna = [f for f in files if f.endswith('.fna.gz') and not f.startswith('cds')]
        gff = [f for f in files if f.endswith('.gff.gz')]
        if not fna or not gff:
            return
        fna.sort(); gff.sort()
        fasta_path = os.path.join(org_dir, fna[0])
        gff_path = os.path.join(org_dir, gff[0])

        # parse fasta
        try:
            fasta_file = gzip.open(fasta_path, 'rt')
        except Exception as e:
            print(f"[{organism}] cannot open fasta: {e}")
            return

        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            if 'mitochondrion' in seq_record.description:
                continue
            if "chromosome" not in seq_record.description:
                continue
            ch_id = seq_record.id
            if str(ch_id) in existing_completed:
                # skip pre-existing rows
                continue
            seq = str(seq_record.seq)
            if len(seq)/1e6 < 0.05:
                continue

            # masking
            exon_mask_seq, Intron_P, Intron_C, Exon_P, Exon_C, CDS_P, CDS_C = exon_and_intron_masking(seq, gff_path, ch_id)
            # skip if no CDS in either strand
            if CDS_P == 0 or CDS_C == 0:
                continue

            # nucleotide analysis
            df_ann = nucleotide_analysis_for_exon(exon_mask_seq)

            # build the row dict with all 134 columns
            row = build_row_from_counts(organism, ch_id, seq, exon_mask_seq, df_ann, (Intron_P, Intron_C, Exon_P, Exon_C, CDS_P, CDS_C))
            # push to writer queue
            queue.put(row)
            # free memory
            try:
                del exon_mask_seq, df_ann, row
            except Exception:
                pass
            print(f"[{organism}] queued {ch_id}")

        fasta_file.close()
    except Exception as e:
        print(f"Error processing organism {organism}: {e}")
    return

# ----------------- main controller -----------------
def gene_element_annotation_in_plus_and_minis(all_genomes_folder_path, results_file, num_workers=None):
    t0 = time.time()
    if num_workers is None:
        num_workers = multiprocessing.cpu_count()
    else:
        num_workers = int(num_workers)

    # prepare existing completed Chr_IDs if results_file exists (resume)
    existing_completed = set()
    if os.path.exists(results_file):
        try:
            df_done = pd.read_csv(results_file, usecols=["Chr_ID"])
            existing_completed = set(df_done["Chr_ID"].astype(str).tolist())
            print(f"Found {len(existing_completed)} already processed chromosomes — will resume and skip them.")
        except Exception:
            existing_completed = set()

    # create manager queue and start writer process
    manager = Manager()
    queue = manager.Queue(maxsize=1000)  # buffer up to 1000 rows
    writer = Process(target=writer_proc, args=(queue, results_file, COLUMN_ORDER, existing_completed), daemon=True)
    writer.start()

    # start worker pool
    organisms = [o for o in os.listdir(all_genomes_folder_path) if not (o.endswith(".json") or o.endswith(".jsonl"))]
    print(f"Starting processing: {len(organisms)} organisms, {num_workers} workers.")
    with ProcessPoolExecutor(max_workers=num_workers) as exe:
        futures = []
        for org in organisms:
            futures.append(exe.submit(process_one_organism, org, all_genomes_folder_path, queue, existing_completed))
        # wait for completion
        for fut in as_completed(futures):
            try:
                fut.result()
            except Exception as e:
                print("Worker exception:", e)

    # all workers done -> send sentinel to writer and wait
    queue.put(None)
    writer.join(timeout=3600)
    elapsed = (time.time() - t0) / 60.0
    print("Finished in {:.2f} minutes".format(elapsed))
    return

# ----------------- run -----------------
if __name__ == "__main__":
    gene_element_annotation_in_plus_and_minis(ONE_TRABITE, RESULTS_FILE, NUM_WORKERS)
    print("Results saved to:", RESULTS_FILE)
