
all_genomes_folder_path = r"/media/su/Bioinf/Genomes_2025_3/Genomes_2025_3_7/Fungi_annotated_compelet_chromosome_genome_gff_cds_2025/ncbi_dataset/data/" # pc 5118
print('inalysis_started')
# listing file
import os
import gzip
import pandas as pd
import numpy as np
import time
start_time = time.time()
from datetime import datetime
current_time = datetime.now()
# Format it as a string for naming files (e.g., "2024-12-15_14-30-45")
start_run = current_time.strftime("%Y-%m-%d_%H-%M-%S")
print(start_run)
#Intron
A2T_list =                      []
G2C_list =                      []
P2C_number_list =               [[],[],[],[],[]]
GC_countent_list_W =            []
GC_countent_list_P =            [[],[],[],[],[]]
GC_countent_list_C =            [[],[],[],[],[]]
len_seq_list =                  []
plus_and_minis_strand_portion = [[],[],[],[],[]]
organims_name =                 []
A2T_anot_list =                 [[],[],[],[],[]]
G2C_anot_list =                 [[],[],[],[],[]]
P_element_number_list =         [[],[],[],[],[]]
C_element_number_list =         [[],[],[],[],[]]

P_A_usage_list =                [[],[],[],[],[]]
P_T_usage_list =                [[],[],[],[],[]]
P_C_usage_list =                [[],[],[],[],[]]
P_G_usage_list =                [[],[],[],[],[]]

C_A_usage_list =                [[],[],[],[],[]]
C_T_usage_list =                [[],[],[],[],[]]
C_C_usage_list =                [[],[],[],[],[]]
C_G_usage_list =                [[],[],[],[],[]]
smaller_portion_2_Bigger =      [[],[],[],[],[]]
dif2sum =                       [[],[],[],[],[]]
P_Portion_list =                [[],[],[],[],[]]
C_portion_list =                [[],[],[],[],[]]
P2C_list =                      [[],[],[],[],[]]
P2C_list_W =                    [[],[],[],[],[]]

A2T_P_list =                    [[],[],[],[],[]]
G2C_P_list =                    [[],[],[],[],[]]
A2T_C_list =                    [[],[],[],[],[]]
G2C_C_list =                    [[],[],[],[],[]]

chr_id_list =                   []
All_annotation =                []
Int_annotaion =                 []
UTR_annotation =                []
CDS_annotation =                []
Shared_annotation =             []
Un_annotation =                 []
Sum_E =                         []

Gene_P_counter_list  =               []
Gene_C_counter_list  =               []
Exon_P_counter_list  =               []
Exon_C_counter_list  =               []
CDS_P_counter_list  =               []
CDS_C_counter_list  =               []


analysed_chr_counter = 0
lablename = 'refrence'
organ_group = "fungis"

organism_list = os.listdir(all_genomes_folder_path)
#for organism in organisms:
organism_counter = 0
#print(organism_list)


for organism in organism_list: # 1 time use not in loop

    if organism.endswith(".json") or organism.endswith(".jsonl"):
        continue
    organism_counter+=1
    #print(organism,' analysis')
    files = os.listdir(all_genomes_folder_path+'/'+organism)
    #print(files)

    fna = [i for i in files if i.endswith('.fna.gz') and not i.startswith('cds')]
    #out = [j for j in files if j.endswith('.out')]
    gffs = [k for k in files if k.endswith('.gff.gz')]
    csv = [p for p in files if p.endswith('.csv')]
    fna.sort()
    #out.sort()
    gffs.sort()
    from all_functions import  exon_and_intron_masking
    from all_functions import opening_fasta_file
    from all_functions import nucleotide_analysis_for_exon
    from Bio import SeqIO


    try:
        fasta_file = gzip.open(all_genomes_folder_path+'/'+organism+'/'+fna[0],'rt')
    except IOError or EOFError:
        continue
    try:
        analysed_chr_counter = 0
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            if 'mitochondrion' in seq_record.description:
                continue
                #print(seq_record.description)
            analysed_chr_counter += 1
            sequence = str(seq_record.seq)
            l  = len(sequence)/1000000
            if l < 0.05:
                continue
            if "chromosome" not in seq_record.description:
                continue
            chr_id = seq_record.id
            chr_number = (seq_record.description.split())[0]
            #print('sequence id :',ch_id,' len chr :',l,'Mb')
            #print('opening_fasta_file finished')

            gff_path  = all_genomes_folder_path+'/'+organism+'/'+gffs[0]


            chr_id = chr_id.split('.')[0]
            sequence = sequence.lower()  # atcg all
            sequence_plus = sequence
            sequence_plus = list(sequence_plus) #atcg all

            Inreon_P_counter = 0
            Inreon_C_counter = 0
            UTR_P_counter = 0
            UTR_C_counter = 0
            CDS_P_counter = 0
            CDS_C_counter = 0
            #lableing genes***************************************************************************
            if gff_path.endswith('.gz'):
                gff = gzip.open(gff_path)
            else:
                gff = open(gff_path)
            df = pd.read_csv(gff,sep='\t',comment='#')
            df.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
            counter2 = 0
            ch_id_split = chr_id.split('.')[0]

            filtered_df = df[(df['seqid'].str.startswith(ch_id_split)) & (df['type'].str.startswith('gene'))]

            for annotation in filtered_df.to_numpy():

                counter2 += 1
                start = int(annotation[3])
                end = int(annotation[4])
                if end - start == 0:
                    #print('len was 0')
                    continue
                if annotation[6] == '+':
                    Inreon_P_counter += 1
                    gene_sequence = sequence_plus[(start - 1):end]
                    gene_sequence = ''.join(gene_sequence)
                    gene_sequence = gene_sequence.replace('a','Q')
                    gene_sequence = gene_sequence.replace('t','W')
                    gene_sequence = gene_sequence.replace('c','E')
                    gene_sequence = gene_sequence.replace('g','R')
                    gene_sequence = gene_sequence.replace('q','A')
                    gene_sequence = gene_sequence.replace('w','T')
                    gene_sequence = gene_sequence.replace('e','C')
                    gene_sequence = gene_sequence.replace('r','G')
                    sequence_plus[(start - 1):end] = gene_sequence

                if annotation[6] == '-':
                    Inreon_C_counter += 1
                    gene_sequence = sequence_plus[(start - 1):end]
                    gene_sequence = ''.join(gene_sequence)
                    gene_sequence = gene_sequence.replace('a','q')
                    gene_sequence = gene_sequence.replace('t','w')
                    gene_sequence = gene_sequence.replace('c','e')
                    gene_sequence = gene_sequence.replace('g','r')
                    gene_sequence = gene_sequence.replace('Q','A')
                    gene_sequence = gene_sequence.replace('W','T')
                    gene_sequence = gene_sequence.replace('E','C')
                    gene_sequence = gene_sequence.replace('R','G')
                    sequence_plus[(start - 1):end] = gene_sequence

            gff.close()
            gff.close()
            #print('gene pass',len(sequence_plus))

            #lableing exons****************************************************************************

            if gff_path.endswith('.gz'):
                gff = gzip.open(gff_path)
            else:
                gff = open(gff_path)

            df = pd.read_csv(gff,sep='\t',comment='#')
            df.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
            counter2 = 0
            ch_id_split = chr_id.split('.')[0]
            filtered_df = df[(df['seqid'].str.startswith(ch_id_split)) & (df['type'].str.startswith('exon'))]
            for annotation in filtered_df.to_numpy():
                counter2 += 1
                start = int(annotation[3])
                end = int(annotation[4])
                if end - start == 0:
                    #print('len was 0')
                    continue

                if annotation[6] == '+':
                    UTR_P_counter += 1
                    exon_sequence = sequence_plus[(start - 1):end]
                    exon_sequence = ''.join(exon_sequence)
                    exon_sequence = exon_sequence.replace('Q','S')
                    exon_sequence = exon_sequence.replace('W','D')
                    exon_sequence = exon_sequence.replace('E','F')
                    exon_sequence = exon_sequence.replace('R','H')
                    exon_sequence = exon_sequence.replace('q','A')
                    exon_sequence = exon_sequence.replace('w','T')
                    exon_sequence = exon_sequence.replace('e','C')
                    exon_sequence = exon_sequence.replace('r','G')
                    sequence_plus[(start - 1):end] = exon_sequence

                if annotation[6] == '-':
                    UTR_C_counter += 1
                    exon_sequence = sequence_plus[(start - 1):end]
                    exon_sequence = ''.join(exon_sequence)
                    exon_sequence = exon_sequence.replace('q','s')
                    exon_sequence = exon_sequence.replace('w','d')
                    exon_sequence = exon_sequence.replace('e','f')
                    exon_sequence = exon_sequence.replace('r','h')
                    exon_sequence = exon_sequence.replace('Q','A')
                    exon_sequence = exon_sequence.replace('W','T')
                    exon_sequence = exon_sequence.replace('E','C')
                    exon_sequence = exon_sequence.replace('R','G')
                    sequence_plus[(start - 1):end] = exon_sequence
            gff.close()
            gff.close()
            #print('exon pass',len(sequence_plus))


            #CDS ******************* CDS ********************* CDS
            if gff_path.endswith('.gz'):
                gff = gzip.open(gff_path)
            else:
                gff = open(gff_path)

            df = pd.read_csv(gff,sep='\t',comment='#')
            df.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
            counter2 = 0
            ch_id_split = chr_id.split('.')[0]
            filtered_df = df[(df['seqid'].str.startswith(ch_id_split)) & (df['type'].str.startswith('CDS'))]
            for annotation in filtered_df.to_numpy():
                counter2 += 1
                start = int(annotation[3])
                end = int(annotation[4])
                if end - start == 0:
                    #print('len was 0')
                    continue

                if annotation[6] == '+':
                    CDS_P_counter+=1
                    CDS_sequence = sequence_plus[(start - 1):end]
                    CDS_sequence = ''.join(CDS_sequence)
                    CDS_sequence = CDS_sequence.replace('S','Z')
                    CDS_sequence = CDS_sequence.replace('D','X')
                    CDS_sequence = CDS_sequence.replace('F','V')
                    CDS_sequence = CDS_sequence.replace('H','B')
                    CDS_sequence = CDS_sequence.replace('z','A')
                    CDS_sequence = CDS_sequence.replace('x','T')
                    CDS_sequence = CDS_sequence.replace('v','C')
                    CDS_sequence = CDS_sequence.replace('b','G')
                    sequence_plus[(start - 1):end] = CDS_sequence

                if annotation[6] == '-':
                    CDS_C_counter +=1
                    CDS_sequence = sequence_plus[(start - 1):end]
                    CDS_sequence = ''.join(CDS_sequence)
                    CDS_sequence = CDS_sequence.replace('s','z')
                    CDS_sequence = CDS_sequence.replace('d','x')
                    CDS_sequence = CDS_sequence.replace('f','v')
                    CDS_sequence = CDS_sequence.replace('h','b')
                    CDS_sequence = CDS_sequence.replace('Z','A')
                    CDS_sequence = CDS_sequence.replace('X','T')
                    CDS_sequence = CDS_sequence.replace('V','C')
                    CDS_sequence = CDS_sequence.replace('B','G')
                    sequence_plus[(start - 1):end] = CDS_sequence
            gff.close()
            gff.close()
            #print('CDS pass',len(sequence_plus))

            sequence_plus = ''.join(sequence_plus)
            #print('return joining',len(sequence_plus))
            exon_mask_seq,Gene_P_counter,Gene_C_counter,Exon_P_counter,Exon_C_counter =sequence_plus ,Inreon_P_counter,Inreon_C_counter,UTR_P_counter,UTR_C_counter

            #*******************************************************************************************
            if CDS_P_counter == 0 or CDS_C_counter == 0:
                continue


            Gene_P_counter_list.append(Gene_P_counter)
            Gene_C_counter_list.append(Gene_C_counter)
            Exon_P_counter_list.append(Exon_P_counter)
            Exon_C_counter_list.append(Exon_C_counter)
            CDS_P_counter_list.append(CDS_P_counter)
            CDS_C_counter_list.append(CDS_C_counter)

            A_W =  sequence.count('A')+sequence.count("a")
            T_W =  sequence.count('T')+sequence.count("t")
            C_W =  sequence.count('C')+sequence.count("c")
            G_W =  sequence.count("G")+sequence.count("g")
            N_W =  sequence.count('N')+sequence.count("n")
            Length = len(sequence) - N_W
            organims_name.append(organism)
            GC_countent_list_W.append(round((G_W+C_W)/(A_W+T_W+C_W)*100,4))
            chr_id_list.append(chr_id)
            len_seq_list.append(round(Length/1000000,6))
            A2T_list.append(round((A_W)/(T_W),4))
            G2C_list.append(round((G_W)/(C_W),4))

            #print('exon_and_intron_masking finished')

            L = len(exon_mask_seq)

            #******************************************************************************************************

            sequenced_length = len(exon_mask_seq)

            entron_P = round((exon_mask_seq.count('Q') + exon_mask_seq.count('W') + exon_mask_seq.count('E') + exon_mask_seq.count('R'))/sequenced_length*100,2)
            exon_p =   round((exon_mask_seq.count('S') + exon_mask_seq.count('D') + exon_mask_seq.count('F') + exon_mask_seq.count('H'))/sequenced_length*100,2)
            CDS_P =    round((exon_mask_seq.count('Z') + exon_mask_seq.count('X') + exon_mask_seq.count('V') + exon_mask_seq.count('B'))/sequenced_length*100,2)
            entron_C = round((exon_mask_seq.count('q') + exon_mask_seq.count('w') + exon_mask_seq.count('e') + exon_mask_seq.count('r'))/sequenced_length*100,2)
            exon_C =   round((exon_mask_seq.count('s') + exon_mask_seq.count('d') + exon_mask_seq.count('f') + exon_mask_seq.count('h'))/sequenced_length*100,2)
            CDS_C =    round((exon_mask_seq.count('z') + exon_mask_seq.count('x') + exon_mask_seq.count('v') + exon_mask_seq.count('b'))/sequenced_length*100,2)
            shered =   round((exon_mask_seq.count('A') + exon_mask_seq.count('T') + exon_mask_seq.count('C') + exon_mask_seq.count('G'))/sequenced_length*100,2)


            #print('plus_information','entron',entron_P,'UTR',exon_p,'CDS',CDS_P)
            #print('comp_information','entron',entron_C,'UTR',exon_C,'CDS',CDS_C,)
            #print('masked_portion',(entron_P+exon_p+CDS_P+entron_C+exon_C+CDS_C),'+ Shered:',shered)


            A_entron_P = exon_mask_seq.count('Q')
            T_entron_P = exon_mask_seq.count('W')
            C_entron_P = exon_mask_seq.count('E')
            G_entron_P = exon_mask_seq.count('R')
            sum_entron_P = A_entron_P + T_entron_P + C_entron_P + G_entron_P
            #Ap_entron_P = round(A_entron_P/sum_entron_P*100,3)
            #Tp_entron_P = round(T_entron_P/sum_entron_P*100,3)
            #Cp_entron_P = round(C_entron_P/sum_entron_P*100,3)
            #Gp_entron_P = round(G_entron_P/sum_entron_P*100,3)

            A_entron_C = exon_mask_seq.count('q')
            T_entron_C = exon_mask_seq.count('w')
            C_entron_C = exon_mask_seq.count('e')
            G_entron_C = exon_mask_seq.count('r')
            sum_entron_C = A_entron_C + T_entron_C + C_entron_C + G_entron_C
            #Ap_entron_C = round(A_entron_C/sum_entron_C*100,3)
            #Tp_entron_C = round(T_entron_C/sum_entron_C*100,3)
            #Cp_entron_C = round(C_entron_C/sum_entron_C*100,3)
            #Gp_entron_C = round(G_entron_C/sum_entron_C*100,3)


            #exon********************
            A_exon_P = exon_mask_seq.count('S')
            T_exon_P = exon_mask_seq.count('D')
            C_exon_P = exon_mask_seq.count('F')
            G_exon_P = exon_mask_seq.count('H')
            sum_exon_P = A_exon_P + T_exon_P + C_exon_P + G_exon_P
            #Ap_exon_P = round(A_exon_P/sum_exon_P*100,3)
            #Tp_exon_P = round(T_exon_P/sum_exon_P*100,3)
            #Cp_exon_P = round(C_exon_P/sum_exon_P*100,3)
            #Gp_exon_P = round(G_exon_P/sum_exon_P*100,3)

            A_exon_C = exon_mask_seq.count('s')
            T_exon_C = exon_mask_seq.count('d')
            C_exon_C = exon_mask_seq.count('f')
            G_exon_C = exon_mask_seq.count('h')
            sum_exon_C = A_exon_C + T_exon_C + C_exon_C + G_exon_C
            #Ap_exon_C = round(A_exon_C/sum_exon_C*100,3)
            #Tp_exon_C = round(T_exon_C/sum_exon_C*100,3)
            #Cp_exon_C = round(C_exon_C/sum_exon_C*100,3)
            #Gp_exon_C = round(G_exon_C/sum_exon_C*100,3)

            #CDS********************
            A_CDS_P = exon_mask_seq.count('Z')
            T_CDS_P = exon_mask_seq.count('X')
            C_CDS_P = exon_mask_seq.count('V')
            G_CDS_P = exon_mask_seq.count('B')
            sum_CDS_P = A_CDS_P + T_CDS_P + C_CDS_P + G_CDS_P
            #Ap_CDS_P = round(A_CDS_P/sum_CDS_P*100,3)
            #Tp_CDS_P = round(T_CDS_P/sum_CDS_P*100,3)
            #Cp_CDS_P = round(C_CDS_P/sum_CDS_P*100,3)
            #Gp_CDS_P = round(G_CDS_P/sum_CDS_P*100,3)

            A_CDS_C = exon_mask_seq.count('z')
            T_CDS_C = exon_mask_seq.count('x')
            C_CDS_C = exon_mask_seq.count('v')
            G_CDS_C = exon_mask_seq.count('b')
            sum_CDS_C = A_CDS_C + T_CDS_C + C_CDS_C + G_CDS_C
            #Ap_CDS_C = round(A_CDS_C/sum_CDS_C*100,3)
            #Tp_CDS_C = round(T_CDS_C/sum_CDS_C*100,3)
            #Cp_CDS_C = round(C_CDS_C/sum_CDS_C*100,3)
            #Gp_CDS_C = round(G_CDS_C/sum_CDS_C*100,3)


            #shared
            A_shared = exon_mask_seq.count('A')
            T_shared = exon_mask_seq.count('T')
            C_shared = exon_mask_seq.count('C')
            G_shared = exon_mask_seq.count('G')
            sum_shared = A_shared + T_shared + C_shared + G_shared
            #Ap_shared = round(A_shared/sum_shared*100,3)
            #Tp_shared = round(T_shared/sum_shared*100,3)
            #Cp_shared = round(C_shared/sum_shared*100,3)
            #Gp_shared = round(G_shared/sum_shared*100,3)


            #unannotated


            A_unannotated = exon_mask_seq.count('a')
            T_unannotated = exon_mask_seq.count('t')
            C_unannotated = exon_mask_seq.count('c')
            G_unannotated = exon_mask_seq.count('g')
            sum_unannotated = A_unannotated + T_unannotated + C_unannotated + G_unannotated
            #Ap_unannotated = round(A_unannotated/sum_unannotated*100,3)
            #Tp_unannotated = round(T_unannotated/sum_unannotated*100,3)
            #Cp_unannotated = round(C_unannotated/sum_unannotated*100,3)
            #Gp_unannotated = round(G_unannotated/sum_unannotated*100,3)





            import pandas as pd
            A = [A_entron_P,A_entron_C,A_exon_P,A_exon_C,A_CDS_P,A_CDS_C,A_shared,A_unannotated]
            T = [T_entron_P,T_entron_C,T_exon_P,T_exon_C,T_CDS_P,T_CDS_C,T_shared,T_unannotated]
            C = [C_entron_P,C_entron_C,C_exon_P,C_exon_C,C_CDS_P,C_CDS_C,C_shared,C_unannotated]
            G = [G_entron_P,G_entron_C,G_exon_P,G_exon_C,G_CDS_P,G_CDS_C,G_shared,G_unannotated]
            dic = {'A':A,'T':T,'C':C,'G':G}
            df = pd.DataFrame(dic)
            df.index = ['+entron','-entron','+UTR','-UTR','+CDS','-CDS','Shared','Unannotated']




            #***********************************************************************************************
            df.loc["Chr_Len"] = {"A":L,'T':L,"C":L,"G":L,"whole":L}
            n_count = exon_mask_seq.count('n')+exon_mask_seq.count('N')
            df.loc["Nn_number"] = {"A":n_count,'T':n_count,"C":n_count,"G":n_count,"whole":n_count}

            df.loc["organism"] = {"A":organism,'T':organism,"C":organism,"G":organism,"whole":organism}
            df.loc["chr_id"] = {"A":chr_id,'T':chr_id,"C":chr_id,"G":chr_id,"whole":chr_id}
            #print(df)


            gene_elements = ['Int','UTR','CDS','Share','Unanot']
            counter = 0
            portions = []
            main_dict={'Name':organims_name,
                       'Chr_ID':chr_id_list,
                       'Length(M)':len_seq_list,
                       'GC_W':GC_countent_list_W,
                       'A2T_W':A2T_list,
                       'G2C_W':G2C_list,}
            for e_number,element in enumerate(gene_elements):
                #print(element)
                if element == "Share" or element == "Unanot":

                    Ap,Tp,Cp,Gp = df.iloc[counter+0]
                    Ap,Tp,Cp,Gp = Ap/2,Tp/2,Cp/2,Gp/2
                    Ac,Tc,Cc,Gc = df.iloc[counter+0]
                    Ac,Tc,Cc,Gc = Ac/2,Tc/2,Cc/2,Gc/2
                    counter += 1
                    #print(element,'+',Ap,Tp,Cp,Gp)
                    #print(element,'-',Ac,Tc,Cc,Gc)
                else:

                    Ap,Tp,Cp,Gp = df.iloc[counter+0]
                    Ac,Tc,Cc,Gc = df.iloc[counter+1]
                    #print(element,'+',Ap,Tp,Cp,Gp)
                    #print(element,'-',Ac,Tc,Cc,Gc)

                    counter += 2

                try :
                    u = 1/Ap + 1/Tp + 1/Cp + 1/Gp
                    u = 1/Ac + 1/Tc + 1/Cc + 1/Gc
                    portions.append(sum([Ap,Tp,Cp,Gp,Ac,Tc,Cc,Gc])/Length*100)

                    #ppppppppppppppppppppppppppppppppp
                    P_portion = round(Ap+Tp+Cp+Gp,4)
                    P_A_usage_list[e_number].append(round(Ap/P_portion*100,4))
                    P_T_usage_list[e_number].append(round(Tp/P_portion*100,4))
                    P_C_usage_list[e_number].append(round(Cp/P_portion*100,4))
                    P_G_usage_list[e_number].append(round(Gp/P_portion*100,4))
                    P_Portion_list[e_number].append(round(P_portion/Length*100,4))
                    A2T_P_list[e_number].append(round(Ap/Tp,4))
                    G2C_P_list[e_number].append(round(Gp/Cp,4))
                    P_element_number_list[e_number].append(CDS_P_counter)
                    GC_countent_list_P[e_number].append(round((Cp+Gp)/P_portion*100,4))

                    #CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                    C_portion = round(Ac+Tc+Cc+Gc,4)
                    C_A_usage_list[e_number].append(round(Ac/C_portion*100,4))
                    C_T_usage_list[e_number].append(round(Tc/C_portion*100,4))
                    C_C_usage_list[e_number].append(round(Cc/C_portion*100,4))
                    C_G_usage_list[e_number].append(round(Gc/C_portion*100,4))
                    C_portion_list[e_number].append(round(P_portion/Length*100,4))
                    A2T_C_list[e_number].append(round(Ac/Tc,4))
                    G2C_C_list[e_number].append(round(Gc/Cc,4))
                    C_element_number_list[e_number].append(CDS_C_counter)
                    GC_countent_list_C[e_number].append(round((Cc+Gc)/C_portion*100,4))

                    #annotat
                    A2T_anot_list[e_number].append(round((Ap+Ac)/(Tp+Tc),4))
                    G2C_anot_list[e_number].append(round((Gp+Gc)/(Cp+Cc),4))


                    P2C_list[e_number].append(round(P_portion/C_portion,4))
                    P2C_number_list[e_number].append(round(CDS_P_counter/CDS_C_counter,4))
                    PandC = [P_portion,C_portion]
                    PandC.sort()
                    smaller_portion_2_Bigger[e_number].append(round(PandC[0]/PandC[1],4))
                    dif2sum[e_number].append(round(abs(PandC[0]-PandC[1])/((PandC[0]+PandC[1])*2)*100,4))
                except ZeroDivisionError:

                    #print(element,'*****',Ap,Tp,Cp,Gp)

                    portions.append(0)
                    #ppppppppppppppppppppppppppppppppp
                    P_portion = round(Ap+Tp+Cp+Gp,4)
                    P_A_usage_list[e_number].append(0)
                    P_T_usage_list[e_number].append(0)
                    P_C_usage_list[e_number].append(0)
                    P_G_usage_list[e_number].append(0)
                    P_Portion_list[e_number].append(0)
                    A2T_P_list[e_number].append(0)
                    G2C_P_list[e_number].append(0)
                    P_element_number_list[e_number].append(0)
                    GC_countent_list_P[e_number].append(0)

                    #CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                    C_portion = round(Ac+Tc+Cc+Gc,4)
                    C_A_usage_list[e_number].append(0)
                    C_T_usage_list[e_number].append(0)
                    C_C_usage_list[e_number].append(0)
                    C_G_usage_list[e_number].append(0)
                    C_portion_list[e_number].append(0)
                    A2T_C_list[e_number].append(0)
                    G2C_C_list[e_number].append(0)
                    C_element_number_list[e_number].append(0)
                    GC_countent_list_C[e_number].append(0)

                    #annotat
                    A2T_anot_list[e_number].append(0)
                    G2C_anot_list[e_number].append(0)

                    P2C_list[e_number].append(0)
                    P2C_number_list[e_number].append(0)
                    PandC = [P_portion,C_portion]
                    PandC.sort()
                    smaller_portion_2_Bigger[e_number].append(0)
                    dif2sum[e_number].append(0)



                vaiavle_dic={
                    f'{element}_GC_P':GC_countent_list_P[e_number],
                    f'{element}_GC_C':GC_countent_list_C[e_number],
                    f'{element}_P%':P_Portion_list[e_number],
                    f'{element}_C%':C_portion_list[e_number],
                    f'{element}_P2C':P2C_list[e_number],
                    f'{element}_S2B':smaller_portion_2_Bigger[e_number],
                    f'{element}_dif2sum':dif2sum[e_number],
                    f'{element}_A2T_ano':A2T_anot_list[e_number],
                    f'{element}_G2C_ano':G2C_anot_list[e_number],
                    f'{element}_P_E_num':P_element_number_list[e_number],
                    f'{element}_C_E_num':C_element_number_list[e_number],
                    f'{element}_P_A_ano':P_A_usage_list[e_number],
                    f'{element}_P_T_ano':P_T_usage_list[e_number],
                    f'{element}_P_C_ano':P_C_usage_list[e_number],
                    f'{element}_P_G_ano':P_G_usage_list[e_number],
                    f'{element}_P_A2T'  :A2T_P_list[e_number],
                    f'{element}_P_G2C'  :G2C_P_list[e_number],
                    f'{element}_C_A_ano':C_A_usage_list[e_number],
                    f'{element}_C_T_ano':C_T_usage_list[e_number],
                    f'{element}_C_C_ano':C_C_usage_list[e_number],
                    f'{element}_C_G_ano':C_G_usage_list[e_number],
                    f'{element}_A2T_C':A2T_C_list[e_number],
                    f'{element}_G2C_C':G2C_C_list[e_number]}

                main_dict = main_dict | vaiavle_dic
            All_annotation.append(round(sum(portions[:4]),4))
            Int_annotaion.append(round(portions[0],4))
            UTR_annotation.append(round(portions[1],4))
            CDS_annotation.append(round(portions[2],4))
            Shared_annotation.append(round(portions[3],4))
            Un_annotation.append(round(portions[4],4))
            Sum_E.append(round(sum(portions),4))
            #print(sum(portions[:4]))
            #main_dict.inser[f'All_anno'] = All_annotation
            try:
                my_table = pd.DataFrame(main_dict)
            except ValueError:
                for i in main_dict:
                    print('hi',i,len(main_dict[i]))
            my_table.insert(4,'Wh_ano',All_annotation)
            my_table.insert(5,'Int_ano',Int_annotaion)
            my_table.insert(6,'UTR_ano',UTR_annotation)
            my_table.insert(7,'CDS_ano',CDS_annotation)
            my_table.insert(8,'share_ano',Shared_annotation)
            my_table.insert(9,'Un_ano',Un_annotation)
            my_table.insert(10,'sum',Sum_E)
            my_table['Gene_P_count'] = Gene_P_counter_list
            my_table['Gene_C_count'] = Gene_C_counter_list
            my_table['Exon_P_count'] = Exon_P_counter_list
            my_table['Exon_C_count'] = Exon_C_counter_list
            my_table['CDS_P_count']  = CDS_P_counter_list
            my_table['CDS_C_count']  = CDS_C_counter_list



            # Print the last row with aligned columns
            column_width = 15# Adjust width as needed
            header = " ".join(col.ljust(column_width) for col in my_table.columns)
            last_row_values = " ".join(str(val).ljust(column_width) for val in my_table.iloc[-1].values)

            # Print the header once
            #if len(my_table) == 1:  # If it's the first row being added
            if analysed_chr_counter == 1:
                print(header)
                print("-" * len(header))  # Add a separator line for clarity

            # Print the last row in aligned format

            print(organism_counter,analysed_chr_counter,last_row_values)
            #print(sum(portions[:4]),portions)

    except EOFError:
        continue
import time
end_time = time.time()

print('procces time',round((end_time-start_time)/60,2),'minit')

current_time = datetime.now()

# Format it as a string for naming files (e.g., "2024-12-15_14-30-45")
timestamp = current_time.strftime("%Y-%m-%d_%H-%M-%S")

results_folder = f"{organ_group}{lablename}_Int_UTR_CDS_{timestamp}"
os.makedirs(results_folder, exist_ok=True)  # Creates the folder if it doesn't exist

# Define the output file path
output_filename = os.path.join(results_folder, f"results_{timestamp}.xlsx")

# Save the DataFrame to the timestamped Excel file
my_table.to_excel(output_filename, index=False)

print(f"Results saved to: {output_filename}")
print('start running',start_run)
print('end running  ',timestamp)
print('procces finished')