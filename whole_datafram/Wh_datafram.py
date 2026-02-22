
'''import pandas as pd


chordata = pd.read_csv(r"C:\Genomics\results\Genes_chardaff_parameters\Chordata.csv")
Viride = pd.read_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Viride_on_saveing_time_resuts.csv')
insect = pd.read_csv(r"C:\Genomics\results\Genes_chardaff_parameters\Arthropoda_500m.csv")
Fungi = pd.read_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Fungisresults_2025-09-17_00-53-54.csv')
Archia = pd.read_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Archaea_2025-11-16_20-58-46.csv')
Bacteria = pd.read_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Bacteria_2025-11-11_17-47-44.csv')
Chloro_plast = pd.read_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Chlroplast_2025-11-18_22-41-13.csv')
Mito = pd.read_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Mito_2025-11-18_19-16-56.csv')
df = pd.concat([chordata, Viride, insect, Fungi, Archia, Bacteria], axis=0, ignore_index=True)




#print(df.head().to_string())


Ap = ((df['Int_P%'] * df['Int_P_A_ano']) +  (df['UTR_P%'] * df['UTR_P_A_ano']) + (df['CDS_P%'] * df['CDS_P_A_ano']) + (df['Share_P%'] * df['Share_P_A_ano']) + (df['Unanot_P%'] * df['Unanot_P_A_ano']))/100
Ac = ((df['Int_C%'] * df['Int_C_A_ano']) +  (df['UTR_C%'] * df['UTR_C_A_ano']) + (df['CDS_C%'] * df['CDS_C_A_ano']) + (df['Share_C%'] * df['Share_C_A_ano']) + (df['Unanot_C%'] * df['Unanot_C_A_ano']))/100
Tp = ((df['Int_P%'] * df['Int_P_T_ano']) +  (df['UTR_P%'] * df['UTR_P_T_ano']) + (df['CDS_P%'] * df['CDS_P_T_ano']) + (df['Share_P%'] * df['Share_P_T_ano']) + (df['Unanot_P%'] * df['Unanot_P_T_ano']))/100
Tc = ((df['Int_C%'] * df['Int_C_T_ano']) +  (df['UTR_C%'] * df['UTR_C_T_ano']) + (df['CDS_C%'] * df['CDS_C_T_ano']) + (df['Share_C%'] * df['Share_C_T_ano']) + (df['Unanot_C%'] * df['Unanot_C_T_ano']))/100

Ap = ( (df['CDS_P%'] * df['CDS_P_A_ano']))/100
Ac = ( (df['CDS_C%'] * df['CDS_C_A_ano']))/100
Tp = ( (df['CDS_P%'] * df['CDS_P_T_ano']))/100
Tc = ( (df['CDS_C%'] * df['CDS_C_T_ano']))/100






b = df[['Int_P%' , 'Int_P_A_ano','UTR_P%', 'UTR_P_A_ano','CDS_P%','CDS_P_A_ano','Share_P%','Share_P_A_ano','Unanot_P%','Unanot_P_A_ano']]
#print('sum_annotation',df['Int_P%']+df['UTR_P%']+df['CDS_P%']+df['Share_P%']+df['Unanot_P%']+(df['Int_C%']+df['UTR_C%']+df['CDS_C%']+df['Share_C%']+df['Unanot_C%']))

#  inputs
A2T = df["CDS_P_A2T"]
P2C = df["CDS_P2C"]
L = df["Length(M)"]



My_formula = ((df["CDS_P2C"] * df["CDS_P_A2T"]) + 1) / (df["CDS_P2C"] + df["CDS_P_A2T"])
A_W_predict =   My_formula/(My_formula+1) * L
T_W_predict =  (1- My_formula/(My_formula+1)) * L
A2T_W_predected = A_W_predict/T_W_predict



df['A2T_coding_information'] = (Ap+Tc)/(Tp+Ac)
My_formula_average_A2T = ((df["CDS_P2C"] *df['A2T_coding_information']+ 1 )) / (df["CDS_P2C"] + df['A2T_coding_information'])

df['full_calculation'] =  (Ap+Ac)/(Tp+Tc)
#print(df['A/T'].to_string())

import matplotlib.pyplot as plt
#plt.scatter(df['A/T'], df['A/T'])
#plt.scatter(df['new_A2T'], df['new_A2T'])


plt.hist(df['CDS_A2T_ano'],bins=400,color = 'black')#csv calculted
plt.hist(My_formula,bins=400,color = 'green') #my M formula
#plt.hist(My_formula_average_A2T,bins=400,color = 'yellow')# my furmola by full A2T
#plt.hist(df['full_calculation'],bins=400,color = 'red') #calculatio real csv amount by calculation
plt.show()'''



# ****************************************************************************************************8
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

'''Bacteria = pd.read_csv(r"C:\Genomics\results\Genes_chardaff_parameters\Bacterial_2025-09-23_20-16-24.csv")
Bacteria['Group'] = 'Bacteria'
chordata = pd.read_csv(r"C:\Genomics\results\Genes_chardaff_parameters\Chordata.csv")
chordata['Group'] = 'Chordata'
Viride = pd.read_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Viride_on_saveing_time_resuts.csv')
Viride['Group'] = 'Viride'
insect = pd.read_csv(r"C:\Genomics\results\Genes_chardaff_parameters\Arthropoda_500m.csv")
insect['Group'] = 'Insect'
Fungi = pd.read_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Fungisresults_2025-09-17_00-53-54.csv')
Fungi['Group'] = 'Fungi'
Archia = pd.read_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Archaea_2025-11-16_20-58-46.csv')
Archia['Group'] = 'Archia'
Chloro_plast = pd.read_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Chlroplast_2025-11-18_22-41-13.csv')
Chloro_plast['Group'] = 'Chloro_plast'
Mito = pd.read_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Mito_2025-11-18_19-16-56.csv')
Mito['Group'] = 'Mito'
df = pd.concat([chordata,Viride,insect,Fungi,Archia,Bacteria,Chloro_plast,Mito],axis=0,ignore_index=True)'''

df =  pd.read_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\whole_datafram\Whole_datafram_2026.csv',low_memory=False)

df['Genome'] = 'dsDNA'
df = df[df['Wh_ano']>0]

path = r'C:\Genomics\Chargaff\coding\all_in_one_programming\Viruses\Viruses_refrence_Int_UTR_CDS_2026-01-26_16-51-15\Viruses_2026-01-26_16-51-15.csv'
df_v =  pd.read_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\Viruses\Viruses_refrence_Int_UTR_CDS_2026-01-26_16-51-15\Viruses_2026-01-26_16-51-15.csv',low_memory=False)
df_v['Group'] = 'Viruses'
df_v = df_v[df_v['Wh_ano']>0]

df = pd.concat([df,df_v],axis=0,ignore_index=True)

#df = df.drop(columns=[	'Int_P_E_num','Int_C_E_num','UTR_P_E_num','UTR_C_E_num','CDS_P_E_num','CDS_C_E_num','Share_P_E_num','Share_C_E_num','Unanot_P_E_num','Unanot_C_E_num' ])
#print(Bacteria.columns.tolist())

# Find rows where both conditions are True
'''mask = (df['UTR_P%'] == 0) & (df['UTR_C%'] == 0)
mask2 = (df['Share_P%'] == 0) & (df['Share_C%'] == 0)
mask3 = (df['Int_P%'] == 0) & (df['Int_C%'] == 0)
mask4 = (df['CDS_P%'] == 0) & (df['CDS_C%'] == 0)
# Set 'UTR_P2C' to 0 only for those rows
df.loc[mask, 'UTR_P2C'] = 0
df.loc[mask2, 'Share_P2C'] = 0
df.loc[mask3, 'Int_P2C'] = 0
df.loc[mask4, 'CDS_P2C'] = 0


maxraj_portion =   df['Int_ano'] +  df['UTR_ano']  + df['CDS_ano'] + df['share_ano']

A2T_P = (df["Int_P_A2T"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_P_A2T"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_P_A2T"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_P_A2T"]*(df['share_ano']/maxraj_portion))
A2T_C = (df["Int_C_A2T"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_C_A2T"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_C_A2T"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_C_A2T"]*(df['share_ano']/maxraj_portion))

G2C_P = (df["Int_P_G2C"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_P_G2C"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_P_G2C"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_P_G2C"]*(df['share_ano']/maxraj_portion))
G2C_C = (df["Int_C_G2C"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_C_G2C"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_C_G2C"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_C_G2C"]*(df['share_ano']/maxraj_portion))


P2C_W = (df["Int_P2C"]   * (df['Int_ano']/maxraj_portion)) + (df["UTR_P2C"]  *(df['UTR_ano']/maxraj_portion)) +(df["CDS_P2C"]  *(df['CDS_ano']/maxraj_portion)) + (df["Share_P2C"]  *(df['share_ano']/maxraj_portion))

#  (df['CDS_P_A2T'] * (df['CDS_P%']/(df['CDS_P%']+df['CDS_C%']))) + (df['CDS_C_A2T'] * (df['CDS_C%']/(df['CDS_P%']+df['CDS_C%'])))

A2T_annotation = (df["Int_A2T_ano"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_A2T_ano"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_A2T_ano"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_A2T_ano"]*(df['share_ano']/maxraj_portion))
G2C_annotation = (df["Int_G2C_ano"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_G2C_ano"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_G2C_ano"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_G2C_ano"]*(df['share_ano']/maxraj_portion))

Ap_annotation = (df["Int_P_A_ano"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_P_A_ano"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_P_A_ano"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_P_A_ano"]*(df['share_ano']/maxraj_portion))
Tp_annotation = (df["Int_P_T_ano"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_P_T_ano"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_P_T_ano"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_P_T_ano"]*(df['share_ano']/maxraj_portion))
Cp_annotation = (df["Int_P_C_ano"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_P_C_ano"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_P_C_ano"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_P_C_ano"]*(df['share_ano']/maxraj_portion))
Gp_annotation = (df["Int_P_G_ano"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_P_G_ano"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_P_G_ano"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_P_G_ano"]*(df['share_ano']/maxraj_portion))

Ac_annotation = (df["Int_C_A_ano"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_C_A_ano"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_C_A_ano"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_C_A_ano"]*(df['share_ano']/maxraj_portion))
Tc_annotation = (df["Int_C_T_ano"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_C_T_ano"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_C_T_ano"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_C_T_ano"]*(df['share_ano']/maxraj_portion))
Cc_annotation = (df["Int_C_C_ano"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_C_C_ano"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_C_C_ano"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_C_C_ano"]*(df['share_ano']/maxraj_portion))
Gc_annotation = (df["Int_C_G_ano"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_C_G_ano"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_C_G_ano"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_C_G_ano"]*(df['share_ano']/maxraj_portion))


df["annot_P_A"] = Ap_annotation
df["annot_P_T"] = Tp_annotation
df["annot_P_C"] = Cp_annotation
df["annot_P_G"] = Gp_annotation

df["annot_C_A"] = Ac_annotation
df["annot_C_T"] = Tc_annotation
df["annot_C_C"] = Cc_annotation
df["annot_C_G"] = Gc_annotation

df['A2T_P'] = A2T_P
df['A2T_C'] = A2T_C

df['G2C_P'] = G2C_P
df['G2C_C'] = G2C_C


df['A2T_annot'] = A2T_annotation
df['G2C_annot'] = G2C_annotation

df['P2C_W'] = P2C_W #annotated


df['ATp_Skew'] = (Ap_annotation - Tc_annotation)/(Ap_annotation + Tp_annotation)
df['GCp_Skew'] = (Gp_annotation - Cp_annotation)/(Gp_annotation + Cp_annotation)
df['ATc_Skew'] = (Ac_annotation - Tc_annotation)/(Ap_annotation + Tc_annotation)
df['GCc_Skew'] = (Gc_annotation - Cc_annotation)/(Gc_annotation + Cc_annotation)
df['ATpc_Skew'] = ((Ap_annotation + Ac_annotation) - (Tp_annotation + Tc_annotation))/((Ap_annotation + Ac_annotation) + (Tp_annotation + Tc_annotation))
df['GCpc_Skew'] = ((Gp_annotation + Gc_annotation) - (Cp_annotation + Cc_annotation))/((Gp_annotation + Gc_annotation) + (Cp_annotation + Cc_annotation))



df['A2T_P-A2T_C'] = df['A2T_P'] - df['A2T_C']
df['G2C_P-G2C_C'] = df['G2C_P'] - df['G2C_C']'''


print(df.isnull().sum())



#df.to_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\whole_datafram\Whole_datafram_2026.csv',index=False)
#df = pd.read_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\whole_datafram\Whole_datafram_2026.csv')


df.to_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\whole_datafram\Whole_datafram_2026_plus_viruse.csv',index=False)
df = pd.read_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\whole_datafram\Whole_datafram_2026_plus_viruse.csv')


#df =df[(df['Chr_ID'].str.startswith('NC'))&(df['Wh_ano']>10)]
df = df.sort_values('A2T_P')

A2T_P = df['A2T_P']
x = np.arange(len(A2T_P))

plt.scatter(x,A2T_P)
plt.show()





