import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from multiprocessing import Pool, cpu_count
import matplotlib
from functools import partial

def violon(plot_df):
    okabe_ito = ['springgreen', '#00FF00', '#00FFFF', '#FFFF00', 'red',
                 'crimson', '#00FF00', '#00FFFF', '#FFFF00', 'red']

    # Set matplotlib to use a serif font for Nature-style plotting
    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['font.serif'] = ['Times New Roman'] + matplotlib.rcParams['font.serif']
    matplotlib.rcParams['font.size'] = 10
    matplotlib.rcParams['axes.linewidth'] = 0.8
    matplotlib.rcParams['lines.linewidth'] = 0.8
    matplotlib.rcParams['legend.fontsize'] = 8
    matplotlib.rcParams['xtick.labelsize'] = 8
    matplotlib.rcParams['ytick.labelsize'] = 8

    # Set style for Nature journal
    sns.set(style="white", font="Times New Roman", font_scale=1.0)

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(8.5, 5), dpi=600)

    # Create violin plot
    sns.violinplot(data=plot_df, palette=okabe_ito, inner="box", cut=0, ax=ax, linewidth=0.6)

    # Customize plot for Nature style
    ax.set_title("Distribution of Sequence Portions and Nucleotide Percentages(Mitochonderia)", fontsize=12, pad=15, weight='bold')
    ax.set_xlabel("", fontsize=10)
    ax.set_ylabel("Percentage (%)", fontsize=10)

    # Customize x-axis labels
    ax.set_xticklabels([
        'P', 'A+', 'T+', 'C+', 'G+',
        'M', 'A-', 'T-', 'C-', 'G-'
    ], rotation=45, ha='right',fontsize=8,color='black',weight='bold')

    # Add mean and std above each violin
    for i, column in enumerate(plot_df.columns):
        mean_val = plot_df[column].mean()
        std_val = plot_df[column].std()
        max_val = plot_df[column].max()
        text = f'μ={mean_val:.1f}\nσ={std_val:.1f}'
        above = 5
        side = 0
        if i == 0  :
            above = -7
            side = 0.8
        if i == 5 :
            above = -7
            side = 0.4
        ax.text(i+side, max_val + above, text, ha='center', va='bottom', fontsize=8)

    # Adjust spines and grid
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.yaxis.grid(True, linestyle='--', linewidth=0.5, color='gray')
    ax.xaxis.grid(False)

    # Adjust layout to prevent label cutoff
    plt.tight_layout()

    # Save and show plot
    plt.savefig("Violons/Mitochonderia_violin_plot_nature_style.png", dpi=900, bbox_inches='tight', format='png')
    plt.show()








if os.path.exists("/media/su/EC3A11743A113D50/genome_analysis/fungi/Violons/Mitochonderia_violon_10_multi.csv"):
    df = pd.read_csv("/media/su/EC3A11743A113D50/genome_analysis/fungi/Violons/Mitochonderia_violon_10_multi.csv")
    print(df.shape)
    violon(df)
else:
    print('it is as first analysis')



    # File paths
    sequence_file_path = r"/media/su/Bioinf/Genomes_2025_3/Genomes_2025_3_7/Mitochondrion/ncbi_dataset/data/genome.fna"
    df = pd.read_csv(r"/media/su/Bioinf/Genomes_2025_3/Genomes_2025_3_7/Mitochondrion/ncbi_dataset/data/organelle_gene_report.csv")#), nrows=10000)


    # Define Okabe-Ito color palette (colorblind-friendly)
    okabe_ito = ['#9467bd', '#2ca02c', '#1f77b4', '#ff7f0e', '#d62728',
                 '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

    # Function to process a single record
    def process_record(record, df):
        if not record.id.startswith("NC"):
            return None
        genes_sequence_plus = ''
        genes_sequence_minis = ''
        id = record.id
        sequence = record.seq
        one_mt_df = df[df["accession"] == id]
        plus = 0
        minis = 0
        P_length = 0
        M_length = 0
        for _, row in one_mt_df.iterrows():
            start = row["start"] - 1
            end = row["end"]
            gene_sequence = sequence[start:end]
            length = end - start
            if row["orientation"] == "plus":
                genes_sequence_plus += gene_sequence
                plus += 1
                P_length += length
            else:
                genes_sequence_minis += gene_sequence
                minis += 1
                M_length += length
        if plus == 0 and minis == 0:
            return None
        total_length = len(genes_sequence_plus + genes_sequence_minis)
        P_portion = len(genes_sequence_plus) / total_length * 100 if total_length > 0 else 0
        Ap = genes_sequence_plus.count('A') / P_length * 100 if P_length > 0 else 0
        Tp = genes_sequence_plus.count('T') / P_length * 100 if P_length > 0 else 0
        Cp = genes_sequence_plus.count('C') / P_length * 100 if P_length > 0 else 0
        Gp = genes_sequence_plus.count('G') / P_length * 100 if P_length > 0 else 0
        M_portion = len(genes_sequence_minis) / total_length * 100 if total_length > 0 else 0
        Am = genes_sequence_minis.count('A') / M_length * 100 if M_length > 0 else 0
        Tm = genes_sequence_minis.count('T') / M_length * 100 if M_length > 0 else 0
        Cm = genes_sequence_minis.count('C') / M_length * 100 if M_length > 0 else 0
        Gm = genes_sequence_minis.count('G') / M_length * 100 if M_length > 0 else 0
        return {
            'P_portion': P_portion, 'A_plus': Ap, 'T_plus': Tp, 'C_plus': Cp, 'G_plus': Gp,
            'M_portion': M_portion, 'A_minus': Am, 'T_minus': Tm, 'C_minus': Cm, 'G_minus': Gm
        }

    # Main processing with multiprocessing
    def main():
        # Lists to store results
        results = []

        # Read all records into memory (since SeqIO.parse is not picklable)
        records = list(SeqIO.parse(sequence_file_path, "fasta"))
        print(f"Total records to process: {len(records)}")

        # Create a partial function with df as an argument
        process_func = partial(process_record, df=df)

        # Use multiprocessing
        with Pool(processes=8) as pool:
            results = pool.map(process_func, records)

        # Filter out None results and create lists
        P_portion_s = []
        Aps = []
        Tps = []
        Cps = []
        Gps = []
        M_portion_s = []
        Ams = []
        Tms = []
        Cms = []
        Gms = []

        for result in results:
            if result is not None:
                P_portion_s.append(result['P_portion'])
                Aps.append(result['A_plus'])
                Tps.append(result['T_plus'])
                Cps.append(result['C_plus'])
                Gps.append(result['G_plus'])
                M_portion_s.append(result['M_portion'])
                Ams.append(result['A_minus'])
                Tms.append(result['T_minus'])
                Cms.append(result['C_minus'])
                Gms.append(result['G_minus'])

        # Create a DataFrame for plotting
        data = {
            'P': P_portion_s, 'A_plus': Aps, 'T_plus': Tps, 'C_plus': Cps, 'G_plus': Gps,
            'M': M_portion_s, 'A_minus': Ams, 'T_minus': Tms, 'C_minus': Cms, 'G_minus': Gms
        }
        plot_df = pd.DataFrame(data)
        plot_df.to_csv('Violons/Mitochonderia_violon_10_multi.csv', index=False)
        violon(plot_df)


    if __name__ == "__main__":
        main()