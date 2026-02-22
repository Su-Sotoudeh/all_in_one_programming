import os
import gzip
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
global N

# ----------------- plotting function -----------------
def violon(plot_df):
    print(plot_df.shape,'analysed chr number')
    okabe_ito = ['springgreen', '#00FF00', '#00FFFF', '#FFFF00', 'red',
                 'crimson', '#00FF00', '#00FFFF', '#FFFF00', 'red']

    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['font.size'] = 10
    sns.set(style="white", font="serif", font_scale=1.0)

    fig, ax = plt.subplots(figsize=(8.5, 5), dpi=900)
    sns.violinplot(data=plot_df, palette=okabe_ito, inner="box",
                   cut=0, ax=ax, linewidth=0.6)

    ax.set_title("Distribution of Sequence Portions and Nucleotide Percentages (Fungi)",
                 fontsize=12, pad=15, weight='bold')
    ax.set_ylabel("Percentage (%)", fontsize=10)

    ax.set_xticklabels([
        'P', 'A+', 'T+', 'C+', 'G+',
        'M', 'A-', 'T-', 'C-', 'G-'
    ], rotation=45, ha='right', fontsize=8, color='black', weight='bold')

    for i, column in enumerate(plot_df.columns):
        print(i,column)
        mean_val = plot_df[column].mean()
        std_val = plot_df[column].std()
        max_val = plot_df[column].max()
        text = f'μ={mean_val:.1f}\nσ={std_val:.1f}'
        above = 5
        side = 0
        if i == 0 or i == 5:
            above = -8
            side = 0.4

        ax.text(i+side, max_val + above, text, ha='center', va='bottom', fontsize=8)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.grid(True, linestyle='--', linewidth=0.5, color='gray')

    plt.tight_layout()
    plt.savefig("Fungi_violin_plot_nature_style.png", dpi=600, bbox_inches='tight')
    plt.show()


# ----------------- worker per bacterium -----------------
def process_one_organism(organism, all_genomes_folder_path):

    files = os.listdir(os.path.join(all_genomes_folder_path, organism))
    fna = [f for f in files if (f.endswith('.fna.gz') or f.endswith('.fna')) and not f.startswith('cds')]
    gff = [f for f in files if (f.endswith('.gff.gz') or f.endswith('.gff'))]

    if not fna or not gff:
        return None

    fna_path = os.path.join(all_genomes_folder_path, organism, fna[0])
    gff_path = os.path.join(all_genomes_folder_path, organism, gff[0])

    fasta_file = gzip.open(fna_path, 'rt') if fna_path.endswith('.gz') else open(fna_path, 'rt')
    gff_file = gzip.open(gff_path, 'rt') if gff_path.endswith('.gz') else open(gff_path, 'r')

    try:
        gff_df = pd.read_csv(gff_file, sep='\t', comment='#',
                             names=['seqid','source','type','start','end','score','strand','phase','attributes'])
    except Exception:
        return None

    results = []
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        if 'plasmid' in seq_record.description or 'mitochondrion' in seq_record.description:
            continue
        ch_id = seq_record.id
        if ch_id.startswith('NW') or "chromosome" not in seq_record.description:
            continue

        sequence = str(seq_record.seq).upper()
        filtered_df = gff_df[gff_df['seqid'] == ch_id]
        filtered_df = filtered_df[filtered_df['type'] == 'gene']

        genes_sequence_plus, genes_sequence_minis = '', ''
        plus, minis, P_length, M_length = 0, 0, 0, 0

        for _, row in filtered_df.iterrows():
            start, end = int(row["start"]) - 1, int(row["end"])
            gene_sequence = sequence[start:end]
            length = end - start
            if row["strand"] == "+":
                genes_sequence_plus += gene_sequence
                plus += 1
                P_length += length
            elif row["strand"] == "-":
                genes_sequence_minis += gene_sequence
                minis += 1
                M_length += length

        if plus + minis == 0:
            continue

        total_length = len(genes_sequence_plus + genes_sequence_minis)
        data = {
            'P': len(genes_sequence_plus) / total_length * 100 if total_length > 0 else 0,
            'A_plus': genes_sequence_plus.count('A') / P_length * 100 if P_length > 0 else 0,
            'T_plus': genes_sequence_plus.count('T') / P_length * 100 if P_length > 0 else 0,
            'C_plus': genes_sequence_plus.count('C') / P_length * 100 if P_length > 0 else 0,
            'G_plus': genes_sequence_plus.count('G') / P_length * 100 if P_length > 0 else 0,
            'M': len(genes_sequence_minis) / total_length * 100 if total_length > 0 else 0,
            'A_minus': genes_sequence_minis.count('A') / M_length * 100 if M_length > 0 else 0,
            'T_minus': genes_sequence_minis.count('T') / M_length * 100 if M_length > 0 else 0,
            'C_minus': genes_sequence_minis.count('C') / M_length * 100 if M_length > 0 else 0,
            'G_minus': genes_sequence_minis.count('G') / M_length * 100 if M_length > 0 else 0,
            'Fungi': organism
        }
        results.append(data)

    fasta_file.close()
    gff_file.close()
    return pd.DataFrame(results)


# ----------------- controller -----------------
def gene_element_annotation_in_plus_and_minis(all_genomes_folder_path):
    refrence = '/media/su/Bioinf/Genomes_2025_3/Genomes_2025_3_7/Fungi_annotated_compelet_chromosome_genome_gff_cds_2025/Fungis_refrence_annot_chr_2025_10_17.tsv'
    Z = pd.read_csv(refrence, sep='\t')
    Z = Z["Assembly Accession"].tolist()

    organism_list = [o for o in os.listdir(all_genomes_folder_path)]# if o in Z]
    print(len(organism_list),'organisms shuld be annalysed number')
    max_workers = 8#min(8, multiprocessing.cpu_count())
    results = []

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_org = {executor.submit(process_one_organism, org, all_genomes_folder_path): org
                         for org in organism_list}
        for future in as_completed(future_to_org):
            try:
                res = future.result()
                if isinstance(res, pd.DataFrame):
                    results.append(res)
            except Exception as e:
                print(f"Error processing {future_to_org[future]}: {e}")

    if results:
        final_df = pd.concat(results, ignore_index=True)
        final_df.to_csv("Fungi_violin_plot_nature_style.csv", index=False)
        return final_df
    else:
        return None


# ----------------- main -----------------
if os.path.exists("Fungi_violin_plot_nature_style.csv"):
    print("CSV found, skipping analysis...")
    df = pd.read_csv("Fungi_violin_plot_nature_style.csv")
else:
    print("Running first-time analysis...")
    base_folder = r"/media/su/Bioinf/Genomes_2025_3/Genomes_2025_3_7/Fungi_annotated_compelet_chromosome_genome_gff_cds_2025/ncbi_dataset/data"
    df = gene_element_annotation_in_plus_and_minis(base_folder)

if df is not None:
    violon(df.drop(columns=["Fungi"]))


print(N)