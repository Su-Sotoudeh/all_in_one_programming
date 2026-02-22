import gzip
import time
start_time = time.time()
from datetime import datetime

current_time = datetime.now()

# Format it as a string for naming files (e.g., "2024-12-15_14-30-45")
start_run = current_time.strftime("%Y-%m-%d_%H-%M-%S")
print(start_run)

import matplotlib
from Bio import SeqIO
import pandas as pd
import numpy as np

from  func_for_fungi import informational_composition_of_organism
from func_for_fungi import seq_Nucleotide_composition
from func_for_fungi import my_function_for_sorting_two_list
from func_for_fungi import equality_usage_of_strand
import matplotlib.pyplot as plt
from func_for_fungi import whole_chromosome_nucleotide_chargaff_test

from func_for_fungi import  whole_annoteted_chromosome_nucleotide_chargaff_test

#whole
A_compusition_list = []
T_compusition_list = []
C_compusition_list = []
G_compusition_list = []
A_T_compusition_list  = [] # A minis T
G_C_compusition_list = [] # G minis C
A2T_list = []
C2G_list = []
P2C_number_list =[]
GC_countent_list = []
GC_countent_list_w = []
GC_countent_list_P = []
GC_countent_list_C = []
len_seq_list = []
portion_of_maping_list = []
plus_and_minis_strand_portion = []
organims_name = []
P_n_usage_list = []
C_n_usage_list = []
A2T_anot_list = []
C2G_anot_list = []
P_element_number_list = []
C_element_number_list = []

P_A_usage_list = []
P_T_usage_list = []
P_C_usage_list = []
P_G_usage_list = []

C_A_usage_list = []
C_T_usage_list = []
C_C_usage_list = []
C_G_usage_list = []

Ap2Tc_list = []
Tp2Ac_list = []
Cp2Gc_list = []
Gp2Cc_list = []
Ac2Tp_list = []
Tc2Ap_list = []
Cc2Gp_list = []
Gc2Cp_list = []

smaller_portion_2_Bigger = []
dif2sum =[]

P_Portion_list = []
C_portion_list = []
Non_portion_list = []
sum_portions_list = []

A2T_P_list = []
C2G_P_list = []
A2T_C_list = []
C2G_C_list = []

N_A_usage_list = []
N_T_usage_list = []
N_C_usage_list = []
N_G_usage_list = []
A2T_n_list = []
C2G_n_list = []





element = 'CDS'

analysed_chr_counter = 0
lablename = 'refrence'


import re
import os
#last_foder_path = r'M:\Genomics\Genomes\fungi2\ncbi_dataset\data'
last_foder_path = r"/media/su/Bioinf/Genomes_2025_3/Genomes_2025_3_7/Fungi_annotated_compelet_chromosome_genome_gff_cds_2025/ncbi_dataset/data/"
organism_list = os.listdir(last_foder_path)
#refrence_list = os.listdir(r"/media/su/Bioinf/Genomes_2025_3/Genomes_2025_3_7/Fungi_annotated_compelet_chromosome_genome_gff_cds_2025/)
import pandas as pd
refrence_list = pd.read_csv(r"/media/su/Bioinf/Genomes_2025_3/Genomes_2025_3_7/Fungi_annotated_compelet_chromosome_genome_gff_cds_2025/fungi_ref_annotated_2025_10_14.tsv",sep='\t')
refrence_list = list(refrence_list['Assembly Accession'])
print(refrence_list)
organism_counter = 0
for organism_directery in organism_list:
    file_list = []
    '''if organism_counter == 2:
        break'''
    #print(organism_directery,refrence_list)
    if organism_directery not in refrence_list:
       continue
    organism_files = os.listdir(fr'{last_foder_path}{organism_directery}')
    organism_files = os.listdir(fr'{last_foder_path}{organism_directery}')
    for file in organism_files:
        files = fr'{last_foder_path}/{file}'
        file_list.append(files)
    if len(file_list) < 2:
        continue



    chr_counter = 0
    #print('****',organism_files)
    for file in organism_files:
        if file.startswith('cds_'):
            continue
        #print(file)
        if file.endswith('.gff.gz'):
            print(file)
            gff_file_path = fr'{last_foder_path}/{organism_directery}/{file}'
            gff_file = gzip.open(gff_file_path, 'rt')
            #gff_file = open(gff_file_path)
            file_content = gff_file.read()
            map_lines = file_content.split('\n')
            if len(map_lines) < 50 :
                print('gff Erroooooooooooooooooor')
                continue
            gff_file.close()

        if file.endswith('.fna.gz') and 'cds_from_genomic' not in file and 'rna.fna.gz' not in file:
            sequence_file_path = file
            print(file)

    organism_counter +=1

    try:
        #zipfile = open(fr'M:\Genomics\Genomes\fungi2\ncbi_dataset\data\{organism_directery}\{sequence_file_path}')
        #zipfile = open(fr'{last_foder_path}{organism_directery}/{sequence_file_path}')
        zipfile = gzip.open(fr'{last_foder_path}{organism_directery}/{sequence_file_path}','rt')

    #print(f'new organism {organism_counter} {organism_directery}')
        for chr in  SeqIO.parse(zipfile,'fasta'):
            #print(chr.id)
            if len(chr.seq) < 50000 :
                continue
            sequence = str(chr.seq)
            len_sequene = (len(sequence))
            chr_number = (chr.description.split())[0]
            annoteted_genome_sequence,P_element_number,C_element_number = equality_usage_of_strand(sequence,gff_file_path,element,chr_number)
            #print(P_element_number,C_element_number)
            if P_element_number == 0 or C_element_number == 0:
                continue
            analysed_chr_counter+=1
            if 'chromosome' in chr.description :
                organims_name.append(f'{organism_directery}(Chr)')
            else:
                organims_name.append(f'{organism_directery}')


            chr_counter +=1





            Ap,Tp,Cp,Gp,Ac,Tc,Cc,Gc,An,Tn,Cn,Gn,N_sum = whole_annoteted_chromosome_nucleotide_chargaff_test(annoteted_genome_sequence)
            #print(organims_name[-1],len_sequene,len(map_lines),Ap,Tp,Cp,Gp,Ac,Tc,Cc,Gc,An,Tn,Cn,Gn,N_sum)
            #******[Ap,Tp,Cp,Gp,Ac,Tc,Cc,Gc,An,Tn,Cn,Gn,N_sum]*************************************************


            #ppppppppppppppppppppppppppppppppp
            P_portion = round(Ap+Tp+Cp+Gp,4)
            P_A_usage_list.append(round(Ap/P_portion*100,4))
            P_T_usage_list.append(round(Tp/P_portion*100,4))
            P_C_usage_list.append(round(Cp/P_portion*100,4))
            P_G_usage_list.append(round(Gp/P_portion*100,4))

            P_Portion_list.append(P_portion)
            Ap2Tc_list.append(round(Ap/Tc,4))
            Tp2Ac_list.append(round(Tp/Ac,4))
            Cp2Gc_list.append(round(Cp/Gc,4))
            Gp2Cc_list.append(round(Gp/Cc,4))
            A2T_P_list.append(round(Ap/Tp,4))
            C2G_P_list.append(round(Cp/Gp,4))
            P_element_number_list.append(P_element_number)



            #CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            C_portion = round(Ac+Tc+Cc+Gc,4)
            C_A_usage_list.append(round(Ac/C_portion*100,4))
            C_T_usage_list.append(round(Tc/C_portion*100,4))
            C_C_usage_list.append(round(Cc/C_portion*100,4))
            C_G_usage_list.append(round(Gc/C_portion*100,4))

            C_portion_list.append(C_portion)
            Ac2Tp_list.append(round(Ac/Tp,4))
            Tc2Ap_list.append(round(Tc/Ap,4))
            Cc2Gp_list.append(round(Cc/Gp,4))
            Gc2Cp_list.append(round(Gc/Cp,4))
            A2T_C_list.append(round(Ac/Tc,4))
            C2G_C_list.append(round(Cc/Gc,4))
            C_element_number_list.append(C_element_number)



            #NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
            N_portion = round(An+Tn+Cn+Gn,4)
            N_A_usage_list.append(round(An/N_portion*100,4))
            N_T_usage_list.append(round(Tn/N_portion*100,4))
            N_C_usage_list.append(round(Cn/N_portion*100,4))
            N_G_usage_list.append(round(Gn/N_portion*100,4))

            Non_portion_list.append(N_portion)
            A2T_n_list.append(round(Ac/Tc,4))
            C2G_n_list.append(round(Cc/Gc,4))



            A2T_anot_list.append(round((Ap+Ac)/(Tp+Tc),4))
            C2G_anot_list.append(round((Cp+Cc)/(Gp+Gc),4))
            sum_portions_list.append(round(N_sum,4))
            GC_countent_list_w.append(round((Cp+Gp+Cc+Gc+Cn+Gn),4))
            GC_countent_list_P.append(round((Cp+Gp)/P_portion*100,4))
            GC_countent_list_C.append(round((Cc+Gc)/C_portion*100,4))
            A2T_list.append(round((Ap+Ac+An)/(Tp+Tc+Tn),4))
            C2G_list.append(round((Cp+Cc+Cn)/(Gp+Gc+Gn),4))
            len_seq_list.append(round(len_sequene/1000000,4))


            P2C_number = round(P_element_number/C_element_number,4)
            P2C_number_list.append(P2C_number)
            plus_and_minis_strand_portion.append(round(P_portion/C_portion,4))
            PandC = [P_portion,C_portion]
            PandC.sort()
            smaller_portion_2_Bigger.append(round(PandC[0]/PandC[1],4))
            dif2sum.append(round(abs(PandC[0]-PandC[1])/((PandC[0]+PandC[1])*2)*100,4))


            my_dict = {'Name':organims_name,'Length(M)':len_seq_list,'GC_W':GC_countent_list_w,
                       'GC_P':GC_countent_list_P,'GC_C':GC_countent_list_C,
                       'A2T_W':A2T_list,'C2G_w':C2G_list,
                       'P%':P_Portion_list,'C%':C_portion_list,'Non%':Non_portion_list,'sum':sum_portions_list,
                       'P2C':plus_and_minis_strand_portion,'S2B':smaller_portion_2_Bigger,
                       'dif2sum':dif2sum,
                       'A2T_anotPC_list':A2T_anot_list,'C2G_anotPC_list':C2G_anot_list,
                       'P_E_num':P_element_number_list,
                       'C_E_num':C_element_number_list,
                       'P/C_E_number':P2C_number_list,
                       'P_A_usage_list':P_A_usage_list,'P_T_usage_list':P_T_usage_list,'P_C_usage_list':P_C_usage_list,
                       'P_G_usage_list':P_G_usage_list,
                       'A2T_P':A2T_P_list,'C2G_P':C2G_P_list,
                       'C_A_usage_list':C_A_usage_list,'C_T_usage_list':C_T_usage_list,
                       'C_C_usage_list':C_C_usage_list,'C_G_usage_list':C_G_usage_list,
                       'A2T_C':A2T_C_list,'C2G_C':C2G_C_list,
                       'N_A':N_A_usage_list,'N_T':N_T_usage_list,'N_C':N_C_usage_list,'N_G':N_G_usage_list,
                       'A2T_n':A2T_n_list,'C2G_n':C2G_n_list,
                       'Ap2Tc':Ap2Tc_list,'Tp2Ac':Tp2Ac_list,'Cp2Gc':Cp2Gc_list,'Gp2Cc':Gp2Cc_list,
                       'Ac2Tp':Ac2Tp_list,'Tc2Ap':Tc2Ap_list,'Cc2Gp':Cc2Gp_list,'Gc2Cp':Gc2Cp_list}




            try:
                my_table = pd.DataFrame(my_dict)
            except ValueError:
                for i in my_dict:
                    print('hi',i,len(my_dict[i]))

            # Print the last row with aligned columns
            column_width = 22 # Adjust width as needed
            header = " ".join(col.ljust(column_width) for col in my_table.columns)
            last_row_values = " ".join(str(val).ljust(column_width) for val in my_table.iloc[-1].values)

            # Print the header once
            if len(my_table) == 1:  # If it's the first row being added
                print(header)
                print("-" * len(header))  # Add a separator line for clarity

            # Print the last row in aligned format
            print(analysed_chr_counter,last_row_values)
    except EOFError:
        print('EOFError')
import time
end_time = time.time()

print('procces time',round((end_time-start_time)/60,2),'minit')

current_time = datetime.now()

# Format it as a string for naming files (e.g., "2024-12-15_14-30-45")
timestamp = current_time.strftime("%Y-%m-%d_%H-%M-%S")

results_folder = f"{lablename}_{element}_{timestamp}"
os.makedirs(results_folder, exist_ok=True)  # Creates the folder if it doesn't exist

# Define the output file path
output_filename = os.path.join(results_folder, f"results_{timestamp}.xlsx")

# Save the DataFrame to the timestamped Excel file
my_table.to_excel(output_filename, index=False)

print(f"Results saved to: {output_filename}")

print('start running',start_run)
print('end running  ',timestamp)








