

# listing file
import os
import gzip
import pandas as pd
import numpy as np



one_Trabite = r'C:\Genomics\Chargaff\coding\all_in_one_programming\CSV_output\\' # pc 5118
#Bioinf = '/media/su/Bioinf/output/Test/' # pc 5118
#path  = '/media/hpc-sotoudeh/ADATA_HD720/output/Human/' #hpc
#path = r'D:\PhD\genomes\Mus_19\\' #laptop
path = one_Trabite
organisms = os.listdir(path)
List_Mp  = []
List_Cp  = []
List_Rp  = []
List_RCp = []
List_Mc  = []
List_Cc  = []
List_Rc  = []
List_RCc = []


organims_names = []

organism_counter = 0
for organism in organisms:
    organism_counter += 1
    if organism_counter == 41:
        break
    path = one_Trabite
    print(organism,' analysis')
    path += organism +'\\'
    files = os.listdir(path)

    fna = [i for i in files if i.endswith('.fna')]
    out = [j for j in files if j.endswith('.out')]
    gff = [k for k in files if k.endswith('.gff.gz')]
    csv = [p for p in files if p.endswith('.csv')]
    fna.sort()
    out.sort()
    gff.sort()
    from all_functions import  exon_and_intron_masking
    from all_functions import opening_fasta_file
    from all_functions import nucleotide_analysis_for_exon
    from all_functions import TE_masking_after_exon
    from all_functions import nucleotide_analysis_for_exon_TE
    #print('list_of_chromosome for analysis')
    chr_counter = 0

    for f,o in zip(fna, out):
        chr_counter+=1
        #print('chr:',chr_counter,'  FNA',f,'  OUT:',o,'    GFF:',gff[0])
        if f[:11] != o[:11]:
            print('not equal FNA and OUT')
            print('Eroor')
            import sys
            sys.exit()
    else:
        print('all files are match')
    chr_couner = 0
    for f,o in zip(fna, out):
        chr_couner+=1
        if chr_couner >= 2:
            continue

        ch_id = f[:-4]
        #if (ch_id + '.csv') in csv:
         #   print(ch_id,'befor proccessd and have csv result')
          #  continue
        #print('start new analysis',f,gff[0],o,ch_id)
        #optional

        organims_names.append(organism)
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

        #reverse converting to ATCG and atcg
        '''CDS+UTR+Intron+LINE+SINE+LTR+DNA+
        CDS-UTR-Intron-LINE-SINE-LTR-DNA-Shared_Unannotated_unsequenced'''
        Mp = exon_TE_masked_sequecne.count("***YI")
        Cp  = exon_TE_masked_sequecne.count("YYY*U")
        Rp  = exon_TE_masked_sequecne.count("IY***")
        RCp = exon_TE_masked_sequecne.count("U*YYY")
        Mc  = exon_TE_masked_sequecne.count("+++yi")
        Cc  = exon_TE_masked_sequecne.count("yyy+u")
        Rc  = exon_TE_masked_sequecne.count("iy+++")
        RCc = exon_TE_masked_sequecne.count("u+yyy")



        List_Mp .append(Mp )
        List_Cp .append(Cp )
        List_Rp .append(Rp )
        List_RCp.append(RCp)
        List_Mc .append(Mc )
        List_Cc .append(Cc )
        List_Rc .append(Rc )
        List_RCc.append(RCc)


        print('AAATG',  "+M" ,Mp ,
                        "+C" ,Cp ,
                        "+R" ,Rp ,
                        "+RC",RCp)
        print('aaatg',  "-M" ,Mc ,
                        "-C" ,Cc ,
                        "-R" ,Rc ,
                        "-RC",RCc)
        import pandas as pd
        dictionary = {"+M"  : List_Mp ,
                      "+C"  : List_Cp ,
                      "+R"  : List_Rp ,
                      "+RC" : List_RCp,
                      "-M"  : List_Mc ,
                      "-C"  : List_Cc ,
                      "-R"  : List_Rc ,
                      "-RC" : List_RCc,}

        new_df = pd.DataFrame(dictionary)
        #new_df = new_df[['+M','+RC','-M','-RC']]
        new_df.index = organims_names
        print(new_df.to_string(formatters={col: '{:,.0f}'.format for col in df.columns}))
        #df.to_csv(path+ch_id+'.csv')
        #print('procces finished')
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt

        # Create the heatmap
        plt.figure(figsize=(10, 6))
        sns.heatmap(new_df, annot=True, fmt='d', cmap='coolwarm',yticklabels=new_df.index)  # fmt='d' for integer annotations
        plt.title("Heatmap of Your DataFrame")
        plt.show()
