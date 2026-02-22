import os
import gzip
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
from Bio import SeqIO
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib


# ----------------- plotting function -----------------
def violon(plot_df):
    print(plot_df.shape, 'analysed chr number')
    okabe_ito = ['springgreen', '#00FF00', '#00FFFF', '#FFFF00', 'red',
                 'crimson', '#00FF00', '#00FFFF', '#FFFF00', 'red']

    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['font.size'] = 10
    sns.set(style="white", font="serif", font_scale=1.0)

    fig, ax = plt.subplots(figsize=(8.5, 5), dpi=600)
    sns.violinplot(data=plot_df, palette=okabe_ito, inner="quartile",
                   cut=0, ax=ax, linewidth=0.6)

    ax.set_title("Distribution of Sequence Portions and Nucleotide Percentages (Arthropoda)",
                 fontsize=12, pad=15, weight='bold')
    ax.set_ylabel("Percentage (%)", fontsize=10)

    ax.set_xticklabels([
        'P', 'A+', 'T+', 'C+', 'G+',
        'M', 'A-', 'T-', 'C-', 'G-'
    ], rotation=45, ha='right', fontsize=8, color='black', weight='bold')

    plt.tight_layout()
    plt.savefig("Arthropoda_violin_plot_nature_style.png", dpi=600, bbox_inches='tight')
    plt.show()


# ----------------- worker per chromosome -----------------
def process_one_chromosome(organism, seq_record, gff_df):
    """Analyse one chromosome, return dict with results."""
    if 'plasmid' in seq_record.description or 'mitochondrion' in seq_record.description:
        return None
    ch_id = seq_record.id
    if ch_id.startswith('NW') or "chromosome" not in seq_record.description:
        return None

    sequence = str(seq_record.seq).upper()
    filtered_df = gff_df[gff_df['seqid'] == ch_id]
    filtered_df = filtered_df[filtered_df['type'] == 'gene']

    genes_sequence_plus, genes_sequence_minis = '', ''
    P_length, M_length = 0, 0

    for _, row in filtered_df.iterrows():
        start, end = int(row["start"]) - 1, int(row["end"])
        gene_sequence = sequence[start:end]
        length = end - start
        if row["strand"] == "+":
            genes_sequence_plus += gene_sequence
            P_length += length
        elif row["strand"] == "-":
            genes_sequence_minis += gene_sequence
            M_length += length

    total_length = len(genes_sequence_plus + genes_sequence_minis)
    if total_length == 0:
        return None

    return {
        "Organism": organism,
        "Chr_ID": ch_id,
        "P": len(genes_sequence_plus) / total_length * 100 if total_length > 0 else 0,
        "A_plus": genes_sequence_plus.count('A') / P_length * 100 if P_length > 0 else 0,
        "T_plus": genes_sequence_plus.count('T') / P_length * 100 if P_length > 0 else 0,
        "C_plus": genes_sequence_plus.count('C') / P_length * 100 if P_length > 0 else 0,
        "G_plus": genes_sequence_plus.count('G') / P_length * 100 if P_length > 0 else 0,
        "M": len(genes_sequence_minis) / total_length * 100 if total_length > 0 else 0,
        "A_minus": genes_sequence_minis.count('A') / M_length * 100 if M_length > 0 else 0,
        "T_minus": genes_sequence_minis.count('T') / M_length * 100 if M_length > 0 else 0,
        "C_minus": genes_sequence_minis.count('C') / M_length * 100 if M_length > 0 else 0,
        "G_minus": genes_sequence_minis.count('G') / M_length * 100 if M_length > 0 else 0,
    }


# ----------------- controller -----------------
def gene_element_annotation_in_plus_and_minis(all_genomes_folder_path,
                                              output_csv="Violon_Arthropoda_results.csv",
                                              max_workers=None):

    # Use more workers to hide I/O waits
    if max_workers is None:
        max_workers = multiprocessing.cpu_count() * 2

    # Resume support
    if os.path.exists(output_csv):
        done = pd.read_csv(output_csv, usecols=["Organism", "Chr_ID"])
        done_set = set(zip(done["Arthropoda"], done["Chr_ID"]))
        header_written = True
        print(f"Resuming from {output_csv}, found {len(done_set)} completed chromosomes")
    else:
        done_set = set()
        header_written = False

    organism_list = [o for o in os.listdir(all_genomes_folder_path)
                     if os.path.isdir(os.path.join(all_genomes_folder_path, o))]
    print(len(organism_list), "organisms to be analysed")

    org_total = {}   # total chromosomes per organism
    org_done = defaultdict(int)  # finished chromosomes per organism

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for organism in organism_list:
            org_path = os.path.join(all_genomes_folder_path, organism)
            files = os.listdir(org_path)
            fna = [f for f in files if (f.endswith('.fna.gz') or f.endswith('.fna')) and not f.startswith('cds')]
            gff = [f for f in files if (f.endswith('.gff.gz') or f.endswith('.gff'))]
            if not fna or not gff:
                continue

            fna_path = os.path.join(org_path, fna[0])
            gff_path = os.path.join(org_path, gff[0])
            gff_file = gzip.open(gff_path, 'rt') if gff_path.endswith('.gz') else open(gff_path, 'r')
            gff_df = pd.read_csv(
                gff_file, sep='\t', comment='#',
                names=['seqid','source','type','start','end','score','strand','phase','attributes']
            )
            gff_file.close()

            fasta_file = gzip.open(fna_path, 'rt') if fna_path.endswith('.gz') else open(fna_path, 'rt')
            seq_records = list(SeqIO.parse(fasta_file, "fasta"))
            fasta_file.close()

            org_total[organism] = len(seq_records)

            for seq_record in seq_records:
                key = (organism, seq_record.id)
                if key in done_set:
                    continue
                futures.append(executor.submit(process_one_chromosome, organism, seq_record, gff_df))

        with open(output_csv, "a") as out_f:
            for future in as_completed(futures):
                res = future.result()
                if res:
                    df_row = pd.DataFrame([res])
                    df_row.to_csv(out_f, mode="a", header=not header_written, index=False)
                    header_written = True

                    # progress counter
                    org = res['Arthropoda']
                    org_done[org] += 1
                    if org_done[org] == org_total.get(org, 0):
                        print(f"✅ Finished organism {org} ({len(org_done)} / {len(org_total)})")

    return pd.read_csv(output_csv)


# ----------------- main -----------------
if __name__ == "__main__":
    base_folder = r"/media/su/Bioinf/Genomes_2025_3/Genomes_2025_3_7/Arthropoda_annotated_compelet_chromosome_genome_gff_gff_2025/ncbi_dataset/data"
    csv_file = "Violon_Arthropoda_results.csv"

    df = gene_element_annotation_in_plus_and_minis(base_folder, output_csv=csv_file)

    if df is not None:
        violon(df.drop(columns=["Arthropoda", "Chr_ID"]))
