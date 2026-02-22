import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


genome_path = r"C:\Genomics\Genomes\Bacteria_2025_2_1\bacteria_reference_complet_genomes\ncbi_dataset\data\GCF_000005845.2\GCF_000005845.2_ASM584v2_genomic.fna"
gff_path = r"C:\Genomics\Genomes\Bacteria_2025_2_1\bacteria_reference_complet_genomes\ncbi_dataset\data\GCF_000005845.2\genomic.gff"

from Simple_annotation_lableing import *

sequence = '-AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTG-'


folders = os.listdir(r"C:\Genomics\Genomes\Bacteria_2025_2_1\bacteria_reference_complet_genomes\ncbi_dataset\data\\")
for each_folder in folders:
    if '.json'  in each_folder or '.tsv' in each_folder:
        continue
    path =  os.path.join(r"C:\Genomics\Genomes\Bacteria_2025_2_1\bacteria_reference_complet_genomes\ncbi_dataset\data\\",each_folder)
    genome_path = os.path.join(path,os.listdir(path)[0])
    gff_path = os.path.join(path,os.listdir(path)[1])

    sequence = opening_fasta_file(genome_path)[0]

    annntated_sequence = Plus_minus_lableing(sequence, gff_path)

    sequence = annntated_sequence
    group1 = ['A','T','C','G','a','t','c','g']
    for g in group1:
        for k in group1:
            if k != g:
                sequence = sequence.replace(f'{g}{k}',f'{g}-{k}')


    my_dict = {
        'A_Homopolymer':[],
        'T_Homopolymer':[],
        'C_Homopolymer':[],
        'G_Homopolymer':[],
        'a_Homopolymer':[],
        't_Homopolymer':[],
        'c_Homopolymer':[],
        'g_Homopolymer':[],

    }

    for i in range(1,8):
        for h in group1:
            my_dict[f'{h}_Homopolymer'].append(sequence.count(f'-{h*i}-'))



    df = pd.DataFrame(my_dict,index=np.arange(1,8))
    print(df.to_string())




    plt.style.use('ggplot')
    plt.rcParams['font.size'] = 14
    plt.rcParams['figure.figsize'] = [8,6]
    #plt.plot(my_dict['A_Homopolymer'],color = 'green')
    #plt.plot(my_dict['T_Homopolymer'],color = 'blue')
    plt.plot(my_dict['C_Homopolymer'],color = 'yellow')
    plt.plot(my_dict['G_Homopolymer'],color = 'red')

    #plt.plot(my_dict['a_Homopolymer'],color = 'green')
    #plt.plot(my_dict['t_Homopolymer'],color = 'blue')
    plt.plot(my_dict['c_Homopolymer'],color = 'yellow')
    plt.plot(my_dict['g_Homopolymer'],color = 'red')

    plt.show()