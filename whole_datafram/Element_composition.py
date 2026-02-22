import numpy as np
import pandas as pd
my_datafram = pd.read_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\whole_datafram\Whole_datafram_2026.csv',low_memory=False)

groups = my_datafram['Group'].unique()

for group in groups:
    df = my_datafram[my_datafram['Group'] == group]
    print(group)
    print('organ number',len(df['Name'].unique()))
    print('chr number',len(df['Chr_ID'].unique()))

    L = sum(df['Length(M)'])
    print('whole sequence length is',round(L/1000,3),'Bbp')
    int = sum(df['Int_ano']/100 * df['Length(M)'])/L *100
    utr = sum(df['UTR_ano']/100 * df['Length(M)'])/L *100      
    cds = sum(df['CDS_ano']/100 * df['Length(M)'])/L *100
    shere = sum(df['share_ano']/100 * df['Length(M)'])/L *100
    Un_ano = sum(df['Un_ano']/100 * df['Length(M)'])/L *100

    P_annotaion_portion = sum((df['Int_P%']* df['Length(M)']) + (df['UTR_P%']*df['Length(M)']) + (df['CDS_P%']*df['Length(M)']) + (df['Share_P%']*df['Length(M)']))/100
    C_annotaion_portion = sum((df['Int_C%']* df['Length(M)']) + (df['UTR_C%']* df['Length(M)']) + (df['CDS_C%']* df['Length(M)']) + (df['Share_C%']* df['Length(M)']))/100
    Un_ano_portion = sum(df['Un_ano']/100 * df['Length(M)'])
    print('Unannot_portion',round(Un_ano_portion/1000,3),'Bbp')
    print('+ annotation',round(P_annotaion_portion/1000,3),'Bbp')
    print('- annotation',round(C_annotaion_portion/1000,3),'Bbp')
    print('sum/L',round((Un_ano_portion+P_annotaion_portion+C_annotaion_portion)/L,5))
    print('Int',round(int,3),'%')
    print('UTR',round(utr,3),'%')
    print('CDS',round(cds,3),'%')
    print('Share',round(shere,3),'%')
    print('Un annot',round(Un_ano,3),'%')




