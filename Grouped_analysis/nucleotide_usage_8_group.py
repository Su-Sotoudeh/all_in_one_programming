'''import pandas as pd
Bacteria_file_path = r'C:\\Genomics\results\Genes_chardaff_parameters\Bacterial_2025-09-23_20-16-24.csv'
Fungi_file_path = r'C:\\Genomics\results\Genes_chardaff_parameters\Fungisresults_2025-09-17_00-53-54.csv'
Arthropoda_file_path = r'C:\\Genomics\results\Genes_chardaff_parameters\Arthropoda_500m.csv'
Chordata_file_path = r'C:\\Genomics\results\Genes_chardaff_parameters\Chordata.csv'
Viride_file_path = r'C:\\Genomics\results\Genes_chardaff_parameters\Viride_on_saveing_time_resuts.csv'
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

A = 'CDS_P_A_ano'
T = 'CDS_P_T_ano'
C = 'CDS_P_C_ano'
G = 'CDS_P_G_ano'
A_c = 'CDS_C_A_ano'
A_T = (whole_df[A] - whole_df[T]).tolist()
analysis_df = whole_df
analysis_df['A-T'] = A_T

whole_df.sort_values(['GC_W'],inplace=True)

print(whole_df.head().to_string())
import matplotlib.pyplot as plt
#plt.style.use('ggplot')
#plt.rcParams['font.size'] = 14
#plt.rcParams['axes.labelsize'] = 14
#plt.rcParams['xtick.labelsize'] = 14
#plt.rcParams['ytick.labelsize'] = 14
#plt.rcParams['legend.fontsize'] = 14
#plt.rcParams['figure.figsize'] = (8.0, 4.0)
#plt.rcParams['figure.dpi'] = 300
import numpy as np
x = np.arange(1,len(whole_df[A])+1)
x = x.tolist()
print(x)
plt.scatter(x,whole_df[A],s = 1,color = 'g',label='A')
plt.scatter(x,whole_df[A_c],s = 1,color = 'r',label='A_c')
#plt.scatter(x,whole_df[T],s = 1,color = 'b',label='T')
#plt.scatter(x,whole_df[C],s = 1,color = 'y',label='C')
#plt.scatter(x,whole_df[G],s = 1,color = 'r',label='G')
plt.show()



#test2    ********************************************************************


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# === File paths ===
Bacteria_file_path = r'C:\\Genomics\\results\\Genes_chardaff_parameters\Bacteria_2025-11-11_17-47-44.csv'
Archea_file_path = r'C:\\Genomics\\results\\Genes_chardaff_parameters\Archaea_2025-11-16_20-58-46.csv'
#A = pd.read_excel(Archea_file_path)
#A.to_csv('Archaea_2025-11-16_20-58-46.csv',index=False)
Fungi_file_path = r'C:\\Genomics\\results\\Genes_chardaff_parameters\Fungisresults_2025-09-17_00-53-54.csv'
Arthropoda_file_path = r'C:\\Genomics\\results\\Genes_chardaff_parameters\Arthropoda_500m.csv'
Chordata_file_path = r'C:\\Genomics\\results\\Genes_chardaff_parameters\Chordata.csv'
Viride_file_path = r'C:\\Genomics\\results\\Genes_chardaff_parameters\Viride_on_saveing_time_resuts.csv'
Mito_file_path = r'C:\\Genomics\\results\\Genes_chardaff_parameters\Mito_2025-11-18_19-16-56.csv'
Chlro_file_path = r'C:\\Genomics\\results\\Genes_chardaff_parameters\Chlroplast_2025-11-18_22-41-13.csv'

# === Load CSVs ===
Bacteria_df = pd.read_csv(Bacteria_file_path, index_col=0)
Archia_df = pd.read_csv(Archea_file_path, index_col=0)
Fungi_df = pd.read_csv(Fungi_file_path, index_col=0)
Fungi_df = Fungi_df[Fungi_df['CDS_ano']!=0]
Arthropoda_df = pd.read_csv(Arthropoda_file_path, index_col=0)
Arthropoda_df = Arthropoda_df[Arthropoda_df['CDS_ano']!= 0]
Chordata_df = pd.read_csv(Chordata_file_path, index_col=0)
Chordata_df = Chordata_df[Chordata_df['CDS_ano']!=0]
Viride_df = pd.read_csv(Viride_file_path, index_col=0)
Viride_df = Viride_df[Viride_df['CDS_ano']!=0]
Mito_df = pd.read_csv(Mito_file_path, index_col=0)
Mito_df = Mito_df[Mito_df['CDS_ano']!=0]
Chlro_df = pd.read_csv(Chlro_file_path, index_col=0)
Chlro_df = Chlro_df[Chlro_df['CDS_ano']!=0]

# === Column names ===
Ap = 'CDS_P_A_ano'
Tp = 'CDS_P_T_ano'
Cp = 'CDS_P_C_ano'
Gp = 'CDS_P_G_ano'
Ac = 'CDS_C_A_ano'
Tc = 'CDS_C_T_ano'
Cc = 'CDS_C_C_ano'
Gc = 'CDS_C_G_ano'

# === Sort each dataframe by GC_W ===
for df in [Bacteria_df,Archia_df, Fungi_df, Arthropoda_df, Chordata_df, Viride_df,Mito_df,Chlro_df]:
    df.sort_values(['CDS_GC_P'], inplace=True)

# === Put all into a list for easy looping ===
datasets = [
    ('Bacteria', Bacteria_df),
    ('Archea',Archia_df),
    ('Fungi', Fungi_df),
    ('Arthropoda', Arthropoda_df),
    ('Chordata', Chordata_df),
    ('Viride', Viride_df),
    ('Mito',Mito_df),
    ('Chlro',Chlro_df),
]

# === Create subplots ===
fig, axes = plt.subplots(1, 8, figsize=(18, 2.5), sharey=True)

colors = {'A': 'g', 'T': 'b', 'C': 'y', 'G': 'r'}


# === Plot each dataset ===
for ax, (title, df) in zip(axes, datasets):
    RP = df['CDS_P%']/(df['CDS_P%']+df['CDS_C%'])
    RC = df['CDS_C%']/(df['CDS_P%']+df['CDS_C%'])
    df_A = ((df[Ap]*RP)+(df[Tc]*RC))
    df_T = ((df[Tp]*RP)+(df[Ac]*RC))
    df_C = ((df[Cp]*RP)+(df[Gc]*RC))
    df_G = ((df[Gp]*RP)+(df[Cc]*RC))





    x = df['CDS_GC_P']
    print('Hi',len(x),len(df_A))
    #print(title,df.shape)
    ax.scatter(x, df_A.to_list(), s=1, color=colors['A'], label='A')
    ax.scatter(x, df_T.to_list(), s=1, color=colors['T'], label='T')
    ax.scatter(x, df_C.to_list(), s=1, color=colors['C'], label='C')
    ax.scatter(x, df_G.to_list(), s=1, color=colors['G'], label='G')

    ax.set_title(title)
    ax.set_xlabel('GC %')
    ax.grid(True, alpha=0.3)

axes[0].set_ylabel('Proportion')
axes[0].legend(loc='upper right', markerscale=4, fontsize=8)

plt.tight_layout()
plt.savefig('nucleotide_usage.png',dpi=900)
#plt.show()'''




#test2    ********************************************************************


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# === File paths ===
Bacteria_file_path   = r"C:\Genomics\results\Genes_chardaff_parameters\Bacterial_2025-09-23_20-16-24.csv"
Archea_file_path     = r'C:\\Genomics\results\\Genes_chardaff_parameters\Archaea_2025-11-16_20-58-46.csv'
Fungi_file_path      = r'C:\\Genomics\results\\Genes_chardaff_parameters\Fungisresults_2025-09-17_00-53-54.csv'
Arthropoda_file_path = r'C:\\Genomics\results\\Genes_chardaff_parameters\Arthropoda_500m.csv'
Chordata_file_path   = r'C:\\Genomics\results\\Genes_chardaff_parameters\Chordata.csv'
Viride_file_path     = r'C:\\Genomics\results\\Genes_chardaff_parameters\Viride_on_saveing_time_resuts.csv'
Mito_file_path       = r'C:\\Genomics\results\\Genes_chardaff_parameters\Mito_2025-11-18_19-16-56.csv'
Chlro_file_path      = r'C:\\Genomics\results\\Genes_chardaff_parameters\Chlroplast_2025-11-18_22-41-13.csv'

# === Load CSVs ===
Bacteria_df   = pd.read_csv(Bacteria_file_path, index_col=0)
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
    RP = df['CDS_P%']/(df['CDS_P%']+df['CDS_C%'])
    RC = df['CDS_C%']/(df['CDS_P%']+df['CDS_C%'])
    df_A = ((df[Ap]*RP)+(df[Tc]*RC))
    df_T = ((df[Tp]*RP)+(df[Ac]*RC))
    df_C = ((df[Cp]*RP)+(df[Gc]*RC))
    df_G = ((df[Gp]*RP)+(df[Cc]*RC))

    x = df['CDS_GC_P']

    # Scatter plots
    ax.scatter(x, df_A, s=0.3, color=colors['A'], label='A')
    ax.scatter(x, df_T, s=0.3, color=colors['T'], label='T')
    ax.scatter(x, df_C, s=0.3, color=colors['C'], label='C')
    ax.scatter(x, df_G, s=0.3, color=colors['G'], label='G')

    # Regression lines + stats for each nucleotide
    for y, nuc in zip([df_A, df_T, df_C, df_G], ['A','T','C','G']):
        # Clean x and y before regression
        mask = (~x.isna()) & (~y.isna())
        x_clean, y_clean = x[mask], y[mask]

        if len(x_clean) > 1 and np.std(x_clean) > 0 and np.std(y_clean) > 0:
            slope, intercept, r_value, p_value, std_err = linregress(x_clean, y_clean)
            x_sorted = np.sort(x_clean)
            y_fit = slope * x_sorted + intercept
            ax.plot(x_sorted, y_fit, color='black', linewidth=1)
            ax.text(0.05, 0.95 - 0.05*['A','T','C','G'].index(nuc),
                    f"{nuc}: R²={r_value**2:.2f}",
                    transform=ax.transAxes,
                    fontsize=6, color=colors[nuc], va='top',bbox=dict(boxstyle='round,pad=0.3', fc='beige',alpha=0.6,edgecolor='none'))
        else:
            ax.text(0.05, 0.95 - 0.05*['A','T','C','G'].index(nuc),
                    f"{nuc}: regression not possible",
                    transform=ax.transAxes,
                    fontsize=6, color=colors[nuc], va='top')

    ax.set_title(title, fontsize=9)
    ax.set_xlabel('GC %')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 60)


    
   

axes[0].set_ylabel('Proportion')
axes[0].legend(loc='upper right', markerscale=4, fontsize=8)

plt.tight_layout()
plt.savefig('nucleotide_usage_regression.svg', dpi=900, format='svg',transparent=True)
plt.savefig('nucleotide_usage_regression.png', dpi=3000, format='png',transparent=True)
plt.show()
