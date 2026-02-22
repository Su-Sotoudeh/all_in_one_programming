
import pandas as pd
Bacteria_file_path = r'C:\Genomics\Chargaff\coding\all_in_one_programming\Genes_inf\Bacteria\Bacteria_refrence_Int_UTR_CDS_2025-11-11_17-47-44\Bacteria_2025-11-11_17-47-44.csv'
Fungi_file_path = r'C:\Genomics\results\Genes_chardaff_parameters\Fungisresults_2025-09-17_00-53-54.csv'
Arthropoda_file_path = r'C:\Genomics\results\Genes_chardaff_parameters\Arthropoda_500m.csv'
Chordata_file_path = r'C:\Genomics\results\Genes_chardaff_parameters\Chordata.csv'
Viride_file_path = r'C:\Genomics\results\Genes_chardaff_parameters\Viride_on_saveing_time_resuts.csv'
Bacteria_df = pd.read_csv(Bacteria_file_path,index_col=0)
Fungi_df = pd.read_csv(Fungi_file_path,index_col=0)
Arthropoda_df = pd.read_csv(Arthropoda_file_path,index_col=0)
Chordata_df = pd.read_csv(Chordata_file_path,index_col=0)
Viride_df = pd.read_csv(Viride_file_path,index_col=0)
print(Bacteria_df.shape)
print(Fungi_df.shape)
print(Arthropoda_df.shape)
print(Chordata_df.shape)
print(Viride_df.shape)
whole_df = pd.concat([Bacteria_df,Fungi_df,Arthropoda_df,Chordata_df,Viride_df], axis=0, ignore_index=True)















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
    plt.savefig(r"C:\Genomics\results\Whole_violin_plot_nature_style.png", dpi=900, bbox_inches='tight')
    plt.show()


df1 = pd.read_csv(r"C:\Genomics\Chargaff\coding\all_in_one_programming\Violon_ploting_fileing\CSVs\Bacteria_violin_plot_nature_style.csv")
df2 = pd.read_csv(r"C:\Genomics\Chargaff\coding\all_in_one_programming\Violon_ploting_fileing\CSVs\Archaea_violin_plot_nature_style.csv")
#df3 = pd.read_csv(r"C:\Genomics\Chargaff\coding\all_in_one_programming\Violon_ploting_fileing\CSVs\")
df4 = pd.read_csv(r"C:\Genomics\Chargaff\coding\all_in_one_programming\Violon_ploting_fileing\CSVs\Fungi_violin_plot_nature_style.csv")
df5 = pd.read_csv(r"C:\Genomics\Chargaff\coding\all_in_one_programming\Violon_ploting_fileing\CSVs\chlro_violon_10_multi.csv")
df6 = pd.read_csv(r"C:\Genomics\Chargaff\coding\all_in_one_programming\Violon_ploting_fileing\CSVs\Bacteria_violin_plot_nature_style.csv")
df7 = pd.read_csv(r"C:\Genomics\Chargaff\coding\all_in_one_programming\Violon_ploting_fileing\CSVs\Archaea_violin_plot_nature_style.csv")
df8 = pd.read_csv(r"C:\Genomics\Chargaff\coding\all_in_one_programming\Violon_ploting_fileing\CSVs\Mitochonderia_violon_10_multi.csv")
df = pd.concat([df1, df2, df2, df4,df5.iloc[1:],df6.iloc[1:],df7.iloc[1:],df8.iloc[1:]],ignore_index=True)
df.to_csv('whole_data.csv', index=False)
CDS_P = df['CDS_P%']/(df['CDS_P%'] + df['CDS_C%'])*100
CDS_C = df['CDS_C%']/(df['CDS_P%'] + df['CDS_C%'])*100
df = df[["CDS_P_A_ano",'CDS_P_T_ano','CDS_P_C_ano','CDS_P_G_ano',"CDS_C_A_ano",'CDS_C_T_ano','CDS_C_C_ano','CDS_C_G_ano']]
df.insert(0,'CDS_P%',CDS_P)
df.insert(5,'CDS_C%',CDS_C)

violon(df,'Whole_organims')