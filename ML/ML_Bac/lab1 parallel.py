import gzip
import time
start_time = time.time()
import matplotlib
from Bio import SeqIO

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


element = 'gene'
#A = r'D:\PhD\genomes\bacteria\ncbi_dataset_bac _refrence_annot.zip'
last_shered_foldur = r'C:\Genomics\Genomes\Bacteria_2025_2_1\bacteria_reference_complet_genomes\ncbi_dataset\data'
import os
counter = 0
for file_name in os.listdir(last_shered_foldur):
    if not file_name.startswith('G'):
        continue
    if counter == 10:
        break

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

    #    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    Aw,Tw,Cw,Gw,A2T,C2G = whole_chromosome_nucleotide_chargaff_test(sequence)
    #print('A2T',A2T,'C2G',C2G)
    A2T_list.append(A2T)
    C2G_list.append(C2G)
    #Aw,Tw,Cw,Gw = seq_Nucleotide_composition(sequence)
    #print('Aw',Aw,'Tw',Tw,'Cw',Cw,'Gw',Gw)
    GC_countent_list_w.append(round(Gw+Cw,4))


    '''Plus, Comp, Both = informational_composition_of_organism(sequence, gff_file_path, [element], 50000)

    A, T, C, G = seq_Nucleotide_composition(Both[2])

    A_compusition_list.append(A)
    T_compusition_list.append(T)
    C_compusition_list.append(C)
    G_compusition_list.append(G)
    A_T_compusition_list.append(A-T)
    G_C_compusition_list.append(G-C)
    # not true, GC_countent_list.append(G + C)
    GC_countent_list.append(G + C)
    len_seq_list.append(round(len(sequence)/1000000,3))
    portion_of_maping_list.append(len(Both[2])/len(sequence)*100)'''

    '''if len(Both[2])> len(seq):
        print(Both[2])
        print('Hi')'''
    #2024/10/18
    plus,comp,nonmap,annoteted_genome_sequence,P_n_usage,C_n_usage,P_element_number,C_element_number = equality_usage_of_strand(sequence,gff_file_path,'gene')
    P_Portion_list.append(plus)
    C_portion_list.append(comp)
    Non_portion_list.append(nonmap)
    sum_portions_list.append(round(plus+comp+nonmap,4))

    #print('P_n_usage',P_n_usage)
    #print('C_n_usage',C_n_usage)
    P_A_usage_list.append(P_n_usage[0])
    P_T_usage_list.append(P_n_usage[1])
    P_C_usage_list.append(P_n_usage[2])
    P_G_usage_list.append(P_n_usage[3])
    A2T_P_list.append(round((P_n_usage[0]/P_n_usage[1]),4))
    C2G_P_list.append(round((P_n_usage[2]/P_n_usage[3]),4))

    C_A_usage_list.append(C_n_usage[0])
    C_T_usage_list.append(C_n_usage[1])
    C_C_usage_list.append(C_n_usage[2])
    C_G_usage_list.append(C_n_usage[3])
    A2T_C_list.append(round((C_n_usage[0]/C_n_usage[1]),4))
    C2G_C_list.append(round((C_n_usage[2]/C_n_usage[3]),4))



    P_element_number_list.append(P_element_number)
    C_element_number_list.append(C_element_number)
    P2C_number = round(P_element_number/C_element_number,4)
    P2C_number_list.append(P2C_number)
    plus_and_minis_strand_portion.append(round(plus/comp,4))
    PandC = [plus,comp]
    PandC.sort()
    smaller_portion_2_Bigger.append(round(PandC[0]/PandC[1],4))
    smaller_portion_2_Bigger2.append(round(PandC[0]/(PandC[0]+PandC[1]),4))
    A2T_anot,C2G_anot,P2C_n_L,C2P_n_L,P_C_N_nucleotide_L = whole_annoteted_chromosome_nucleotide_chargaff_test(annoteted_genome_sequence)
    #print('AT2_anot',A2T_anot,'C2G_anot',C2G_anot)
    A2T_anot_list.append(A2T_anot)
    C2G_anot_list.append(C2G_anot)
    #print('Ap2Tc',P2C_n_L[0],'Tp2Ac',P2C_n_L[1],'Cp2Gc',P2C_n_L[2],'Gp2Cc',P2C_n_L[3])
    Ap2Tc_list.append(P2C_n_L[0])
    Tp2Ac_list.append(P2C_n_L[1])
    Cp2Gc_list.append(P2C_n_L[2])
    Gp2Cc_list.append(P2C_n_L[3])


    import numpy as np
    '''x = np.arange(0,len(plus_and_minis_strand_portion),1)
    plt.scatter(len_seq_list,plus_and_minis_strand_portion)
    plt.savefig('CDS_strand usage_equality.svg',transparent=True)
    plt.show()'''


    import pandas as pd

    my_dict = {'Name':organims_name,'Length(M)':len_seq_list,'GC':GC_countent_list_w,'A2T_W':A2T_list,
    'C2G_w':C2G_list,'P%':P_Portion_list,'C%':C_portion_list,'Non%':Non_portion_list,'sum':sum_portions_list
               ,'P2C':plus_and_minis_strand_portion,'S2B':smaller_portion_2_Bigger,
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
               'Ap2Tc':Ap2Tc_list,'Tp2Ac':Tp2Ac_list,'Cp2Gc':Cp2Gc_list,'Gp2Cc':Gp2Cc_list}

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





    #print(my_table.to_string())

my_table.plot(x = 'P2C',y = 'P_E_num',kind = 'scatter')
#plt.show()
import time

my_table.to_excel('my_reult.xlsx')

end_time = time.time()

print('procces time',round(end_time-start_time,2))

from datetime import datetime

# Capture the current date and time
current_time = datetime.now()

# Format it as a string for naming files (e.g., "2024-12-15_14-30-45")
timestamp = current_time.strftime("%Y-%m-%d_%H-%M-%S")


results_folder = f"results_{timestamp}"
os.makedirs(results_folder, exist_ok=True)  # Creates the folder if it doesn't exist

# Define the output file path
output_filename = os.path.join(results_folder, f"results_{timestamp}.xlsx")




# Save the DataFrame to the timestamped Excel file
my_table.to_excel(output_filename, index=False)

# Print confirmation
print(f"Results saved to: {output_filename}")


'''
#new code 2024/06/22
A_T_sort,G_C_sort = my_function_for_sorting_two_list(A_T_compusition_list,G_C_compusition_list)
plt.scatter(GC_countent_list_w,A_compusition_list ,color='green',s = 10)
plt.scatter(GC_countent_list_w,T_compusition_list ,color='blue',s  = 10)
plt.scatter(GC_countent_list_w,C_compusition_list ,color='yellow',s= 10)
plt.scatter(GC_countent_list_w,G_compusition_list ,color='red', s  = 10)
plt.title('Bacteria ' + element)

plt.xlabel('GC %')
plt.ylabel('A,T,C,G %')
plt.savefig('CDS_1.svg',transparent=True)
plt.show()

plt.scatter(GC_countent_list_w,A_T_compusition_list, color = 'green')
plt.scatter(GC_countent_list_w,G_C_compusition_list, color = 'red')
plt.title('Bacteria')
plt.xlabel('GC%')
plt.ylabel('A-T and G-C')
plt.savefig('CDS_2.svg',transparent=True)
plt.show()

plt.plot(A_T_sort, color = 'green')
plt.plot(G_C_sort, color = 'red')
plt.xlabel('max to min diffrence')
plt.ylabel('A-T and G-C')
plt.savefig('CDS_3.svg',transparent=True)
plt.show()

plt.scatter(len_seq_list,GC_countent_list_w, color = 'black')
plt.title('Bacteria')
plt.xlabel('len_sequence')
plt.ylabel('GC %')
plt.savefig('CDS_4.svg',transparent=True)
plt.show()

plt.scatter(portion_of_maping_list,GC_countent_list_w)
plt.title('Bacteria')
plt.xlabel('annotated %')
plt.ylabel('GC %')
plt.savefig('CDS_5.svg',transparent=True)
plt.show()
'''