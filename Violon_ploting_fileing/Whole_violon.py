import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from multiprocessing import Pool, cpu_count
import matplotlib
from functools import partial





def violon(plot_df,viride):
    print(plot_df.shape,'analysed chr number')
    okabe_ito = ['springgreen', '#00FF00', '#00FFFF', '#FFFF00', 'red',
                 'crimson', '#00FF00', '#00FFFF', '#FFFF00', 'red']

    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['font.size'] = 10
    sns.set(style="white", font="serif", font_scale=1.0)

    fig, ax = plt.subplots(figsize=(8.5, 5), dpi=600)
    sns.violinplot(data=plot_df, palette=okabe_ito, inner="box",
                   cut=0, ax=ax, linewidth=0.6)#"quartile"

    ax.set_title("Distribution of Sequence Portions and Nucleotide Percentages (all ogrganism)",
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
        if i == 5 :
            above = -8
            side = 0.4

        ax.text(i+side, max_val + above, text, ha='center', va='bottom', fontsize=8)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.grid(True, linestyle='--', linewidth=0.5, color='gray')

    plt.tight_layout()
    plt.savefig("/media/su/EC3A11743A113D50/genome_analysis/fungi/Whole_violon/Whole_violin_plot_nature_style.png", dpi=900, bbox_inches='tight')
    plt.show()


df1 = pd.read_csv("/media/su/EC3A11743A113D50/genome_analysis/fungi/Viride_refrence_Int_UTR_CDS_CDS_2025-09-19/Viride_on_saveing_time_resuts.csv")
df2 = pd.read_csv("/media/su/EC3A11743A113D50/genome_analysis/fungi/Chordata_refrence_Int_UTR_CDS_2025-09-19_13-22-31/Chordata.csv")
df3 = pd.read_csv("/media/su/EC3A11743A113D50/genome_analysis/fungi/Arthropoda_refrence_Int_UTR_CDS_2025-09-19_13-13-25/Arthropoda_500m.csv")
df4 = pd.read_csv("/media/su/EC3A11743A113D50/genome_analysis/fungi/Fungis_refrence_Int_UTR_CDS_2025-09-17_00-53-54/Fungisresults_2025-09-17_00-53-54.csv")
df5 = pd.read_csv("/media/su/EC3A11743A113D50/genome_analysis/fungi/Chloroplast_violon/chlro_violon_10_multi.csv")
df6 = pd.read_csv("/media/su/EC3A11743A113D50/genome_analysis/fungi/Violons/Bacteria_violin_plot_nature_style.csv")
df7 = pd.read_csv("/media/su/EC3A11743A113D50/genome_analysis/fungi/Violons/Archaea_violin_plot_nature_style.csv")
df8 = pd.read_csv("/media/su/EC3A11743A113D50/genome_analysis/fungi/Violons/Mitochonderia_violon_10_multi.csv")
df = pd.concat([df1, df2, df3, df4,df5.iloc[1:],df6.iloc[1:],df7.iloc[1:],df8.iloc[1:]],ignore_index=True)
df.to_csv('whole_data.csv', index=False)
CDS_P = df['CDS_P%']/(df['CDS_P%'] + df['CDS_C%'])*100
CDS_C = df['CDS_C%']/(df['CDS_P%'] + df['CDS_C%'])*100
df = df[["CDS_P_A_ano",'CDS_P_T_ano','CDS_P_C_ano','CDS_P_G_ano',"CDS_C_A_ano",'CDS_C_T_ano','CDS_C_C_ano','CDS_C_G_ano']]
df.insert(0,'CDS_P%',CDS_P)
df.insert(5,'CDS_C%',CDS_C)

violon(df,'Whole_organims')