one_Trabite = r'C:\Genomics\Chargaff\coding\all_in_one_programming\OUT_CSV_files\\' # pc 5118
#Bioinf = '/media/su/Bioinf/output/Test/' # pc 5118
def all_element_annotation_in_plus_and_minis(genome_folder_path,organism):
    # listing file
    import os
    import gzip
    import pandas as pd
    import numpy as np

    #path  = '/media/hpc-sotoudeh/ADATA_HD720/output/Human/' #hpc
    #path = r'D:\PhD\genomes\Mus_19\\' #laptop
    path = genome_folder_path
    organisms = os.listdir(path)
    print(organisms)

    #for organism in organisms:
    for i in range(1): # 1 time use not in loop
        path = one_Trabite
        print(organism,' analysis')
        path += organism +'/'
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
        print('list_of_chromosome for analysis')
        chr_counter = 0
        for f,o in zip(fna, out):
            chr_counter+=1
            print('chr:',chr_counter,'  FNA',f,'  OUT:',o,'    GFF:',gff[0])
            if f[:11] != o[:11]:
                print('not equal FNA and OUT')
                print('Eroor')
                import sys
                sys.exit()
        else:
            print('all files are match')
        for f,o in zip(fna, out):

            ch_id = f[:-4]
            if (ch_id + '.csv') in csv:
                print(ch_id,'befor proccessd and have csv result')
                continue
            print('start new analysis',f,gff[0],o,ch_id)
            #optional
            seq = opening_fasta_file(path+f)
            print(seq.count('A')+seq.count("a"),seq.count('T')+seq.count("t"),seq.count('C')+seq.count("c"),seq.count("G")+seq.count("g"),seq.count('N')+seq.count("n"),len(seq))
            exon_mask_seq = exon_and_intron_masking(opening_fasta_file(path+f),path+gff[0],ch_id)
            print('exon masking pass')
            #nucleotide_analysis_for_exon(exon_mask_seq)
            exon_TE_masked_sequecne = TE_masking_after_exon(exon_mask_seq,path+o,ch_id)
            print('TE masking pass')
            df = nucleotide_analysis_for_exon_TE(exon_TE_masked_sequecne)
            import pandas as pd
            L = len(exon_TE_masked_sequecne)
            df.loc["Chr_Len"] = {"A":L,'T':L,"C":L,"G":L,"whole":L}
            n_count = exon_TE_masked_sequecne.count('n')+exon_TE_masked_sequecne.count('N')
            df.loc["Nn_number"] = {"A":n_count,'T':n_count,"C":n_count,"G":n_count,"whole":n_count}

            df.loc["organism"] = {"A":organism,'T':organism,"C":organism,"G":organism,"whole":organism}
            df.loc["chr_id"] = {"A":f,'T':f,"C":f,"G":f,"whole":f}


            df.to_csv(path+ch_id+'.csv')
            print('procces finished')
    return organism















