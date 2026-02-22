
import pandas as pd


chordata = pd.read_csv(r"C:\Genomics\results\Genes_chardaff_parameters\Chordata.csv")
Viride = pd.read_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Viride_on_saveing_time_resuts.csv')
insect = pd.read_csv(r"C:\Genomics\results\Genes_chardaff_parameters\Arthropoda_500m.csv")
Fungi = pd.read_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Fungisresults_2025-09-17_00-53-54.csv')
Archia = pd.read_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Archaea_2025-11-16_20-58-46.csv')
Chloro_plast = pd.read_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Chlroplast_2025-11-18_22-41-13.csv')
Mito = pd.read_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Mito_2025-11-18_19-16-56.csv')
df = pd.concat([chordata, Viride, insect, Fungi, Archia], axis=0, ignore_index=True)




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
plt.show()