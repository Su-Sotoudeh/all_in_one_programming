'''import pandas as pd

df = pd.read_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\Viruses\Viruses_refrence_Int_UTR_CDS_2026-01-25_15-56-45\Viruses_2026-01-25_15-56-45.xlsx')
df.to_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\Viruses\Viruses_refrence_Int_UTR_CDS_2026-01-25_15-56-45\Viruses_2026-01-25_15-56-45.csv',index=False)
df.to_excel(r'C:\Genomics\Chargaff\coding\all_in_one_programming\Viruses\Viruses_refrence_Int_UTR_CDS_2026-01-25_15-56-45\Viruses_2026-01-25_15-56-45.xlsx',index=False)

print(df)'''
import os

'''one_Trabite = r"C:\Genomics\Genomes\\unzip_viruse_all_complete\\unzip_viruse_all_complete\\ncbi_dataset\data\\"
organismis = os.listdir(one_Trabite)
print(organismis)
counter = 0
for organism in organismis :
    counter += 1
    if '.jsonl' in organism or 'json' in organism :
        continue
    files = os.listdir(os.path.join(one_Trabite, organism))
    fna = [i for i in files if (i.endswith('.fna') or i.endswith('.fna.gz')) and not i.startswith('cds')]
    gff = [k for k in files if k.endswith('.gff')or k.endswith('.gff.gz')]
    print(counter,len(fna),len(gff))'''



import pandas as pd

df= pd.read_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\Viruses\Viruses_refrence_Int_UTR_CDS_2026-01-26_16-51-15\Viruses_2026-01-26_16-51-15.csv',low_memory=False)

print(df.head().to_string())

df = df[(df['Genome'] == 'ssRNA(-)')&(df['CDS_P_count'] !=0)&(df['CDS_C_count'] != 0)]
print('shape',df.shape)
ids = df['Name'].values
print(ids)

for id in ids:
    one_Trabite = rf"C:\Genomics\Genomes\unzip_viruse_all_complete\unzip_viruse_all_complete\ncbi_dataset\data\{id}"
    organismis = os.listdir(one_Trabite)
    print(organismis)

    
    files = organismis
    fna = [i for i in files if (i.endswith('.fna') or i.endswith('.fna.gz')) and not i.startswith('cds')]
    gff = [k for k in files if k.endswith('.gff')or k.endswith('.gff.gz')]
    gff_file = one_Trabite + '\\' + organismis[1]
    #df_all = pd.read_csv(gff_file)#, sep='\t')#, comment='#')
    file = open(gff_file)
    for line in file:
        print(line)
    #print(df_all.to_string())