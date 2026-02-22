# equality portion between + and - strand

#imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import openpyxl
import gzip
from Bio import SeqIO

# input information
input_organism_name = 'GRCh38.p14'
refseq_Assembly_ID = 'GCF_000001405.40'
input_element_name = ['LINE','SINE','LTR','DNA']#,'RC/Helitron','tRNA','scRNA','Retroposons','RNA','srpRNA'] #study genomic elements
element_name_gff = ['gene']#,'pseudogene','exon','CDS','lnc_RNA','mRNA','enhancer','biological_region','transcript']
chromosomes_names = ['10', '11', '12', '13', '14', '15', '16','17', '18', '19', '20', '21', '22','X','Y']#['10']#['1','2','3','4','5','6','7','8','9' ,'10', '11', '12', '13', '14', '15', '16','17', '18', '19', '20', '21', '22','X','Y']

input_plot_len = 8000
annotation_file_path_out_UCSC = r'D:\PhD\genomes\GRCh38.p14\UCSC\hg38.p14.fa.out'
annotation_file_path_out_refseq = r'D:\PhD\genomes\GRCh38.p14\ncbi_dataset\data\GCF_000001405.40\GCF_000001405.40_GRCh38.p14_rm.out'
fasta_path_refseq = r'D:\PhD\genomes\GRCh38.p14\ncbi_dataset\data\GCF_000001405.40\GCF_000001405.40_GRCh38.p14_genomic.fna'
annotation_file_path_GFF = r'D:\PhD\genomes\GRCh38.p14\ncbi_dataset\data\GCF_000001405.40\genomic.gff'
fasta_path_UCSC = r'D:\PhD\genomes\GRCh38.p14\UCSC\hg38.p14.fa.gz'



from my_functions import  check_and_creat_refseq_file
from functions_for_primats import open_fasta_file_multi_record_UCSC_and_refseq
from functions_for_primats import portion_of_TE_elements_in_geneme_with_UCSC_and_refseq_step2
from functions_for_primats import nucleotide_composition_in_seqenence_simple
from functions_for_primats import nucleotide_composition_in_seqenence_ATCG_2
from functions_for_primats import exon_and_intron_labaling
from functions_for_primats import merge_seq_seq_share_with_1234
from functions_for_primats import analysis_composition_of_labled_sequence
from functions_for_primats import only_ATCG_and_atcg

refseq_records_lists = ['NC_000010.11','NC_000011.10','NC_000012.12','NC_000013.11','NC_000014.9','NC_000015.10','NC_000016.10','NC_000017.11','NC_000018.10','NC_000019.10','NC_000020.11','NC_000021.9','NC_000022.11','NC_000023.11','NC_000024.10']# check_and_creat_refseq_file(refseq_Assembly_ID,input_organism_name,chromosomes_names)


#refseq



#analysis for genes

plus_portions = []
comp_portions = []
duble_portions = []
unmask_portions = []
unsequenced = []


for element in input_element_name:
    for Chr_number,Chr_ID in zip(chromosomes_names,refseq_records_lists):
        Chr = [Chr_number,Chr_ID]
        print(f'start chr{Chr_number,Chr_ID}')
        print('GFF analysis '*10)
        Chr_sequence_refseq = open_fasta_file_multi_record_UCSC_and_refseq(fasta_path_refseq,Chr)


        Plus_strand,Complement_strand,Both = portion_of_TE_elements_in_geneme_with_UCSC_and_refseq_step2(Chr_sequence_refseq,annotation_file_path_out_refseq,Chr,element)
        print('nucleotide_composition_in_seqenence in plus')
        nucleotide_composition_in_seqenence_simple(Plus_strand)
        print('nucleotide_composition_in_seqenence in comp')
        nucleotide_composition_in_seqenence_simple(Complement_strand)
        print('nucleotide_composition_in_seqenence in both')
        nucleotide_composition_in_seqenence_simple(Both)


for Chr_number,Chr_ID in zip(chromosomes_names,refseq_records_lists):
    Chr = [Chr_number,Chr_ID]
    print(f'start chr{Chr_number,Chr_ID}')
    Chr_sequence_refseq = open_fasta_file_multi_record_UCSC_and_refseq(fasta_path_refseq,Chr)
    plus,comp,both = exon_and_intron_labaling(Chr_sequence_refseq,annotation_file_path_GFF,Chr)


















