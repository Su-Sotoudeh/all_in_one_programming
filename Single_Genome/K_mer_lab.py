import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import ConnectionPatch
from matplotlib.colorbar import Colorbar




one_Trabite = r'C:\Genomics\results\OUT_CSV_files\\'  # pc 5118
#Bioinf = '/media/su/Bioinf/output/Test/' # pc 5118
#path  = '/media/hpc-sotoudeh/ADATA_HD720/output/Human/' #hpc
#path = r'D:\PhD\genomes\Mus_19\\' #laptop
path = one_Trabite
organisms = os.listdir(path)



organism_counter = 0


kmers = ['GG','CA','AAT','CCT','ATCG','CAAG','AAATG','GGGGG','AAAAAA','CGGGGA','GCGCGCG','ATATCGC']

kmers = ['AAATG','AAAA','GGGA']

data_dict = {}
for k in kmers:
    for feature in ["+M", "+C", "+R", "+RC"]:
        data_dict[f"{k}_{feature}"] = []



analysed_organisms = []
for organism in organisms:


    #if organism_counter == 41:
    #    break
    path = one_Trabite
    print(organism,' analysis')
    path += organism +'\\'
    files = os.listdir(path)

    fna = [i for i in files if i.endswith('.fna')]
    out = [j for j in files if j.endswith('.out')]
    gff = [k for k in files if k.endswith('.gff.gz')]
    csv = [p for p in files if p.endswith('.csv')]
    print(len(fna),len(out),len(gff),len(csv))
    fna.sort()
    out.sort()
    gff.sort()
    from all_functions import  exon_and_intron_masking
    from all_functions import opening_fasta_file
    from all_functions import nucleotide_analysis_for_exon
    from all_functions import TE_masking_after_exon
    from all_functions import nucleotide_analysis_for_exon_TE
    #print('list_of_chromosome for analysis')

    #print('fna',fna)
    if len(out) == 0:
        print('out',out)

    if len(fna) != len(out) or len(out) == 0:
        continue
    for f,o in zip(fna, out):
        #print('chr:',chr_counter,'  FNA',f,'  OUT:',o,'    GFF:',gff[0])
        if f[:11] != o[:11]:
            print('not equal FNA and OUT')
            print('Eroor')
            import sys
            sys.exit()
    else:
        print('all files are match')
    chr_couner = 0
    organism_counter += 1
    print('start anakysis organism', organism_counter)

    #if organism_counter == 41:
    #    break


    for f,o in zip(fna, out):



        chr_couner+=1

        if chr_couner >= 2:
            con = 0
            break
        print('start chr analysed',chr_couner)
        ch_id = f[:-4]
        #if (ch_id + '.csv') in csv:
        #   print(ch_id,'befor proccessd and have csv result')
        #  continue
        #print('start new analysis',f,gff[0],o,ch_id)
        #optional


        seq = opening_fasta_file(path+f)
        #print('A/T',(seq.count('A')+seq.count("a"))/(seq.count('T')+seq.count("t")),'C/G',(seq.count('C')+seq.count("c"))/(seq.count("G")+seq.count("g")),seq.count('N')+seq.count("n"),len(seq))
        exon_mask_seq = exon_and_intron_masking(opening_fasta_file(path+f),path+gff[0],ch_id)
        #print('exon masking pass')
        #nucleotide_analysis_for_exon(exon_mask_seq)
        exon_TE_masked_sequecne = TE_masking_after_exon(exon_mask_seq,path+o,ch_id)
        #print('TE masking pass')
        list_letters = ['ZXVB','SDFH','QWER','*YUI','OPJK','LM@#','1234',
                        'zxvb','sdfh','qwer','+yui','opjk','lm$%','5678',
                        'ATCG','atcg','~~~~']
        df = nucleotide_analysis_for_exon_TE(exon_TE_masked_sequecne)
        #print(df.to_string(formatters={col: '{:,.2f}'.format for col in df.columns}))

        #return_to_ATCG




        #reverse converting to ATCG and atcg
        '''CDS+UTR+Intron+LINE+SINE+LTR+DNA+
        CDS-UTR-Intron-LINE-SINE-LTR-DNA-Shared_Unannotated_unsequenced'''
        k_mers = ['AAATG','AAAA','GGGA']
        list_k_mers = [["***YI","U*YYY","+++yi","u+yyy"],['****','YYYY','++++','yyyy'],['III*','YUUU','iii+','yuuu']]
        # Suppose we have 2 organisms and 3 k-mers

        for index,l in enumerate(k_mers):
            Mp  = exon_TE_masked_sequecne.count(list_k_mers[index][0])
            RCp = exon_TE_masked_sequecne.count(list_k_mers[index][1])
            Mc  = exon_TE_masked_sequecne.count(list_k_mers[index][2])
            RCc = exon_TE_masked_sequecne.count(list_k_mers[index][3])

            pMp  = round(Mp / (Mp+RCp+Mc+RCc) * 100,2 )
            pRCp = round(RCp/ (Mp+RCp+Mc+RCc) * 100,2 )
            pMc  = round(Mc / (Mp+RCp+Mc+RCc) * 100,2 )
            pRCc = round(RCc/ (Mp+RCp+Mc+RCc) * 100,2 )

            data_dict[f"{l}_+M"].append(pMp)
            data_dict[f"{l}_+RC"].append(pRCp)
            data_dict[f"{l}_-M"].append(pMc)
            data_dict[f"{l}_-RC"].append(pRCc)
            print(l,data_dict)
    analysed_organisms.append(organism)
    print(analysed_organisms)
    print(data_dict)
    df = pd.DataFrame(data_dict, index=analysed_organisms)
    print('df is ...')
    print(df)



    df.index = analysed_organisms
    #df.to_csv('K_mer_heatmap_2_40.csv')
    print(df.to_string(formatters={col: '{:,.0f}'.format for col in df.columns}))

    df.to_csv('140_40.csv')

    #print('procces finished')


    # Create the heatmap
    plt.figure(figsize=(10, 6))
    sns.heatmap(df, annot=True, fmt='2f', cmap='viridis',yticklabels=df.index)  # fmt='d' for integer annotations
    plt.title("Heatmap of Your DataFrame")
    plt.show()
