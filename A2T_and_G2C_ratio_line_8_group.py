
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# === File paths ===
Bacteria_file_path   = r'C:\\Genomics\results\\Genes_chardaff_parameters\Bacterial_2025-09-23_20-16-24.csv'
Archea_file_path     = r'C:\\Genomics\results\\Genes_chardaff_parameters\Archaea_2025-11-16_20-58-46.csv'
Fungi_file_path      = r'C:\\Genomics\results\\Genes_chardaff_parameters\Fungisresults_2025-09-17_00-53-54.csv'
Arthropoda_file_path = r'C:\\Genomics\results\\Genes_chardaff_parameters\Arthropoda_500m.csv'
Chordata_file_path   = r'C:\\Genomics\results\\Genes_chardaff_parameters\Chordata.csv'
Viride_file_path     = r'C:\\Genomics\results\\Genes_chardaff_parameters\Viride_on_saveing_time_resuts.csv'
Mito_file_path       = r'C:\\Genomics\results\\Genes_chardaff_parameters\Mito_2025-11-18_19-16-56.csv'
Chlro_file_path      = r'C:\\Genomics\results\\Genes_chardaff_parameters\Chlroplast_2025-11-18_22-41-13.csv'

# === Load CSVs ===
Bacteria_df   = pd.read_csv(Bacteria_file_path, index_col=0)
print(Bacteria_df.head().to_string())
print(Bacteria_df.shape)
Archia_df     = pd.read_csv(Archea_file_path, index_col=0)
Fungi_df      = pd.read_csv(Fungi_file_path, index_col=0)
Fungi_df      = Fungi_df[Fungi_df['CDS_ano']!=0]
Arthropoda_df = pd.read_csv(Arthropoda_file_path, index_col=0)
Arthropoda_df = Arthropoda_df[Arthropoda_df['CDS_ano']!=0]
Chordata_df   = pd.read_csv(Chordata_file_path, index_col=0)
Chordata_df   = Chordata_df[Chordata_df['CDS_ano']!=0]
Viride_df     = pd.read_csv(Viride_file_path, index_col=0)
Viride_df     = Viride_df[Viride_df['CDS_ano']!=0]
Mito_df       = pd.read_csv(Mito_file_path, index_col=0)
Mito_df       = Mito_df[Mito_df['CDS_ano']!=0]
Chlro_df      = pd.read_csv(Chlro_file_path, index_col=0)
Chlro_df      = Chlro_df[Chlro_df['CDS_ano']!=0]

# === Column names ===
Ap, Tp, Cp, Gp = 'CDS_P_A_ano','CDS_P_T_ano','CDS_P_C_ano','CDS_P_G_ano'
Ac, Tc, Cc, Gc = 'CDS_C_A_ano','CDS_C_T_ano','CDS_C_C_ano','CDS_C_G_ano'

# === Sort each dataframe by GC% ===
for df in [Bacteria_df, Archia_df, Fungi_df, Arthropoda_df, Chordata_df, Viride_df, Mito_df, Chlro_df]:
    df.sort_values(['CDS_GC_P'], inplace=True)

# === Put all into a list for easy looping ===
datasets = [
    ('Bacteria',   Bacteria_df),
    ('Archea',     Archia_df),
    ('Fungi',      Fungi_df),
    ('Arthropoda', Arthropoda_df),
    ('Chordata',   Chordata_df),
    ('Viride',     Viride_df),
    ('Mito',       Mito_df),
    ('Chlor',      Chlro_df),
]

# === Create subplots ===
fig, axes = plt.subplots(1, 8, figsize=(14, 3), sharey=True)

colors = {'A': 'g', 'T': 'b', 'C': 'y', 'G': 'r'}

# === Plot each dataset ===
for ax, (title, df) in zip(axes, datasets):
    RP = df['CDS_P%']/(df['CDS_P%']+df['CDS_C%'])*100
    RC = df['CDS_C%']/(df['CDS_P%']+df['CDS_C%'])*100
    df_A = ((df[Ap]*RP)+(df[Tc]*RC))
    df_T = ((df[Tp]*RP)+(df[Ac]*RC))
    df_C = ((df[Cp]*RP)+(df[Gc]*RC))
    df_G = ((df[Gp]*RP)+(df[Cc]*RC))
    A2T = df_A/df_T
    A2T = A2T/(A2T.max())
    G2C = df_G/df_C
    G2C = G2C/(G2C.max())
    N_df = pd.DataFrame({"A":df_A, "T":df_T, "C":df_C, "G":df_G,'A2T':A2T,'G2C':G2C})
    x = np.arange(len(A2T))
    A =  N_df['A2T'].sort_values()
    G = N_df['G2C'].sort_values()
    print(A)
    # Scatter plots
    ax.scatter( x,A, s=0.3, color=colors['A'], label='A2T')

    ax.scatter( x,G, s=0.3, color=colors['G'], label='G2C')
    #ax.scatter( np.arange(len(RP.tolist())),RP, s=0.3, color=colors['G'], label='RP')

    import numpy as np
    x1 = np.arange(100)
    y1 = np.arange(100)
    #ax.scatter(x1,y1, s=0.3, color='black', label='black line')




    ax.set_title(title, fontsize=9)
    ax.set_xlabel('GC %')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1.5)





axes[0].set_ylabel('Proportion')
axes[0].legend(loc='upper right', markerscale=4, fontsize=8)

plt.tight_layout()
#plt.savefig('nucleotide_usage_line_ratio.svg', dpi=900, format='svg',transparent=True)
#plt.savefig('nucleotide_usage_line_ratio.png', dpi=3000, format='png',transparent=True)
plt.show()