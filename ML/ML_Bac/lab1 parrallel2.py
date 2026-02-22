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

from  func_for_bac import informational_composition_of_organism
from func_for_bac import seq_Nucleotide_composition
from func_for_bac import my_function_for_sorting_two_list
from func_for_bac import equality_usage_of_strand
import matplotlib.pyplot as plt
from func_for_bac import whole_chromosome_nucleotide_chargaff_test

from func_for_bac import  whole_annoteted_chromosome_nucleotide_chargaff_test

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
smaller_portion_2_Bigger2 =[]

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
A = r'D:\PhD\genomes\bacteria\ncbi_dataset_bac _refrence_annot.zip'
last_shered_foldur = r'C:\Genomics\Genomes\Bacteria_2025_2_1\bacteria_reference_complet_genomes\ncbi_dataset\data'
import os
counter = 0
for file_name in os.listdir(last_shered_foldur):
    if not file_name.startswith('G'):
        continue
    '''if counter == 10:
        break'''

    #print(counter)
    bac_foldur = fr'{last_shered_foldur}\{file_name}'
    file_list = []
    for file in os.listdir(bac_foldur):
        files = fr'{bac_foldur}\{file}'
        file_list.append(files)
    #print(file_list)
    if len(file_list) < 2:
        continue
    if file_list[0].endswith('.fna'):
        sequence_path = file_list[0]

    if file_list[1].endswith('.gff'):
        gff_file_path = file_list[1]

    file = open(gff_file_path)
    file_content = file.read()
    map_lines = file_content.split('\n')
    #print(len(map_lines))
    if len(map_lines) < 50 :
        print('gff Erroooooooooooooooooor')
        continue

    for recoed in SeqIO.parse(sequence_path,'fasta') :
        if  not 'plasmid' in recoed.description :
            #print(recoed.description)
            #print('main genoeme')
            #print(recoed.id)
            sequence = recoed.seq
            len_sequene = (len(sequence))
            #print(round(len_sequene/1000000,2),'M')
            break
    if len_sequene < 100000:
        continue
    counter+=1
    organims_name.append(file_name)


    #2024/12/18
    annoteted_genome_sequence,P_element_number,C_element_number = equality_usage_of_strand(sequence,gff_file_path,'CDS')
    Ap,Tp,Cp,Gp,Ac,Tc,Cc,Gc,An,Tn,Cn,Gn,N_sum = whole_annoteted_chromosome_nucleotide_chargaff_test(annoteted_genome_sequence)
    #******[Ap,Tp,Cp,Gp,Ac,Tc,Cc,Gc,An,Tn,Cn,Gn,N_sum]*************************************************


    #ppppppppppppppppppppppppppppppppp
    P_A_usage_list.append(Ap)
    P_T_usage_list.append(Tp)
    P_C_usage_list.append(Cp)
    P_G_usage_list.append(Gp)
    P_portion = round(Ap+Tp+Cp+Gp,4)
    P_Portion_list.append(P_portion)
    Ap2Tc_list.append(round(Ap/Tc,4))
    Tp2Ac_list.append(round(Tp/Ac,4))
    Cp2Gc_list.append(round(Cp/Gc,4))
    Gp2Cc_list.append(round(Gp/Cc,4))
    A2T_P_list.append(round(Ap/Tp,4))
    C2G_P_list.append(round(Cp/Gp,4))
    P_element_number_list.append(P_element_number)



    #CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    C_A_usage_list.append(Ac)
    C_T_usage_list.append(Tc)
    C_C_usage_list.append(Cc)
    C_G_usage_list.append(Gc)
    C_portion = round(Ac+Tc+Cc+Gc,4)
    C_portion_list.append(C_portion)
    Ac2Tp_list.append(round(Ac/Tp,4))
    Tc2Ap_list.append(round(Tc/Ap,4))
    Cc2Gp_list.append(round(Cc/Gp,4))
    Gc2Cp_list.append(round(Gc/Cp,4))
    A2T_C_list.append(round(Ac/Tc,4))
    C2G_C_list.append(round(Cc/Gc,4))
    C_element_number_list.append(C_element_number)



    #NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    N_A_usage_list.append(An)
    N_T_usage_list.append(Tn)
    N_C_usage_list.append(Cn)
    N_G_usage_list.append(Gn)
    N_portion = round(An+Tn+Cn+Gn,4)
    Non_portion_list.append(N_portion)
    A2T_n_list.append(round(Ac/Tc,4))
    C2G_n_list.append(round(Cc/Gc,4))



    A2T_anot_list.append(round((Ap+Ac)/(Tp+Tc),4))
    C2G_anot_list.append(round((Cp+Cc)/(Gp+Gc),4))
    sum_portions_list.append(round(N_sum,4))
    GC_countent_list_w.append(round(Cp+Cc+Gp+Gc,4))
    A2T_list.append(round((Ap+Ac+An)/(Tp+Tc+Tn),4))
    C2G_list.append(round((Cp+Cc+Cn)/(Gp+Gc+Gn),4))
    len_seq_list.append(round(len_sequene/1000000,4))


    P2C_number = round(P_element_number/C_element_number,4)
    P2C_number_list.append(P2C_number)
    plus_and_minis_strand_portion.append(round(P_portion/C_portion,4))
    PandC = [P_portion,C_portion]
    PandC.sort()
    smaller_portion_2_Bigger.append(round(PandC[0]/PandC[1],4))
    smaller_portion_2_Bigger2.append(round(PandC[0]/(PandC[0]+PandC[1]),4))


    my_dict = {'Name':organims_name,'Length(M)':len_seq_list,'GC':GC_countent_list_w,
               'A2T_W':A2T_list,'C2G_w':C2G_list,
               'P%':P_Portion_list,'C%':C_portion_list,'Non%':Non_portion_list,'sum':sum_portions_list,
               'P2C':plus_and_minis_strand_portion,'S2B':smaller_portion_2_Bigger,
               'S2B_2':smaller_portion_2_Bigger2,
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
            print(i,len(my_dict[i]))


    # Print the last row with aligned columns
    column_width = 15  # Adjust width as needed
    header = " ".join(col.ljust(column_width) for col in my_table.columns)
    last_row_values = " ".join(str(val).ljust(column_width) for val in my_table.iloc[-1].values)

    # Print the header once
    if len(my_table) == 1:  # If it's the first row being added
        print(header)
        print("-" * len(header))  # Add a separator line for clarity

    # Print the last row in aligned format
    print(counter,last_row_values)


import time
end_time = time.time()

print('procces time',round((end_time-start_time)/60,2),'minit')

current_time = datetime.now()

# Format it as a string for naming files (e.g., "2024-12-15_14-30-45")
timestamp = current_time.strftime("%Y-%m-%d_%H-%M-%S")

results_folder = f"{element}_{timestamp}"
os.makedirs(results_folder, exist_ok=True)  # Creates the folder if it doesn't exist

# Define the output file path
output_filename = os.path.join(results_folder, f"results_{timestamp}.xlsx")

# Save the DataFrame to the timestamped Excel file
my_table.to_excel(output_filename, index=False)

print(f"Results saved to: {output_filename}")

print('start running',start_run)
print('end running  ',timestamp)

