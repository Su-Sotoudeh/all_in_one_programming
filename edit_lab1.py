import gzip

print("how are you this is my new projects")
import matplotlib
from Bio import SeqIO

from E_coli_functions import CDS_EXON_TE_UNMASKED
from E_coli_functions import CDS_EXON_TE_UNMASKED_plus_minis
from E_coli_functions import plot_CDS_EXON_TE_UNMASKED_plus_minis

from my_functions import  check_and_creat_refseq_file
from functions_for_primats import open_fasta_file_multi_record_UCSC_and_refseq
from functions_for_primats import portion_of_TE_elements_in_geneme_with_UCSC_and_refseq
from functions_for_primats import nucleotide_composition_in_seqenence_simple
from functions_for_primats import nucleotide_composition_in_seqenence_ATCG_2
from functions_for_primats import test
from functions_for_primats import merge_seq_seq_share_with_1234
from functions_for_primats import analysis_composition_of_labled_sequence
from functions_for_primats import only_ATCG_and_atcg







file_locations = r"M:\Genomes\100_genome\\"
import os
genomes_files_list = os.listdir(file_locations)
organism_name1 = 0
organism_name2 = 0
organism_name3 = 0
organim_counter = 0



plus_portions = []
comp_portions = []
duble_portions = []
unmask_portions = []
unsequenced = []
chr_number_list = []
organism_name = []
length_list = []
p2c_list = []



for file in genomes_files_list:

    file_name_spilit = str(file).split(".")
    if file_name_spilit[1] == "fa" and file_name_spilit[2] == "gz":
        fasta_file = file
        organism_name1 = file_name_spilit[0]
        #print(fasta_file)
    if file_name_spilit[1] == "fa" and file_name_spilit[2] == "out":
        out_file = file
        organism_name2 = file_name_spilit[0]
        #print(out_file)
    if file_name_spilit[2] == "gtf":
        gtf_file = file
        organism_name3 = file_name_spilit[0]
        #print(gtf_file)

    if organism_name1 == organism_name2 and organism_name1 == organism_name3:
        organim_counter +=1
        #print(organim_counter)
        #print(organism_name1,organism_name2,organism_name3)
        #print(fasta_file,out_file,gtf_file)
        fasta_path = file_locations + fasta_file
        out_path = file_locations + out_file
        gtf_path = file_locations + gtf_file
        print(fasta_path)
        print(out_path )
        print(gtf_path )
        a = gzip.open(fasta_path,'rt')
        for recoed in SeqIO.parse(a,'fasta') :

            #print(recoed.description)
            if 'mitochondrion'  in recoed.description:
                print('mitochondrion and continue')
                continue
            if 'plasmid' in recoed.description:
                print('plasmid and continue')
                #print(recoed.description)
                continue
            if 'chr' not in recoed.description :
                continue
            if '_NW' in recoed.description :
                continue
            if 'chr_Un' in recoed.description or 'chrUn' in recoed.description \
                    or 'random' in recoed.description or 'chrM' in recoed.description \
                    or 'alt' in recoed.description:
                continue
            Chr_id = recoed.id
            #print(counter,'fungi',Chr_id)
            seq = recoed.seq
            print(organim_counter,recoed.id,len(seq))
            '''gtf_file = gzip.open(gtf_path,'rt')
            for line in gtf_file:
                print(line)
            out_file_open = gzip.open(out_path)
            for line in out_file_open:
                print('out   ',line) '''

            Chr = [str(recoed.id),str(recoed.id)]
            print(f'start chr{(recoed.id + ' ')*50}')
            #print('GFF analysis '*10)
            print('step',1)
            Chr_sequence_refseq = open_fasta_file_multi_record_UCSC_and_refseq(fasta_path,Chr)
            length = len(Chr_sequence_refseq)
            print('step',2)
            plus,comp,both = test(Chr_sequence_refseq,gtf_path,Chr,['transcript','gene'])
            '''print('nucleotide_composition_in_seqenence in plus')
            nucleotide_composition_in_seqenence_simple(plus)
            print('nucleotide_composition_in_seqenence in comp')
            nucleotide_composition_in_seqenence_simple(comp)
            print('nucleotide_composition_in_seqenence in both')
            nucleotide_composition_in_seqenence_simple(both)'''

            #print('GFF finished '+'*'*50)
            #TEs **********************************************
            #print('TEs analysis '+ '*'*50)
            print('step',3)
            Plus_strand,Complement_strand,Both = portion_of_TE_elements_in_geneme_with_UCSC_and_refseq(Chr_sequence_refseq,out_path,Chr)
            '''print('nucleotide_composition_in_seqenence in plus')
            nucleotide_composition_in_seqenence_simple(Plus_strand)
            print('nucleotide_composition_in_seqenence in comp')
            nucleotide_composition_in_seqenence_simple(Complement_strand)
            print('nucleotide_composition_in_seqenence in both')
            nucleotide_composition_in_seqenence_simple(Both)'''


            #print('TEs finished '+'*'*50)
            #nucleotide_composition_in_seqenence_atcg(Plus_strand)
            #nucleotide_composition_in_seqenence_atcg(Complement_strand)
            #nucleotide_composition_in_seqenence_atcg(Both)
            #print('GFF and TEs analysis'+'*' * 50)
            print('step',4)
            merged_sequence_plus = merge_seq_seq_share_with_1234(plus,Plus_strand)
            #print('information about merging of GFF and TEs plus information')
            print('step',5)
            #print('GFF+ and TE+ on plus strand ')
            one_masked_portion_type1_p,one_masked_portion_type2_p,Duble_masked_portion_p,unmasked_region_p,unsequenced_region_p = nucleotide_composition_in_seqenence_ATCG_2(merged_sequence_plus)
            print('step',6)
            plus_sequence = only_ATCG_and_atcg(merged_sequence_plus)
            print('step',7)
            #print('GFF- and TE- on minis strand')
            merged_sequence_comp = merge_seq_seq_share_with_1234(comp,Complement_strand)
            #print('result about merging of GFF and TEs comp information')
            print('step',8)
            one_masked_portion_type1_c,one_masked_portion_type2_c,Duble_masked_portion_c,unmasked_region_c,unsequenced_region_c = nucleotide_composition_in_seqenence_ATCG_2(merged_sequence_comp)
            print('step',9)
            comp_sequence = only_ATCG_and_atcg(merged_sequence_comp)
            print('step',10)
            #print('Both (+-)GFF and (+-)TEs on one strand')
            merged_sequence_both = merge_seq_seq_share_with_1234(both,Both)
            #print('information about merging of GFF and TEs both strand information')
            print('step',11)
            one_masked_portion_type1_B,one_masked_portion_type2_B,Duble_masked_portion_B,unmasked_region_B,unsequenced_region_B = nucleotide_composition_in_seqenence_ATCG_2(merged_sequence_both)
            print('step',12)
            both_sequecne = only_ATCG_and_atcg(merged_sequence_both)
            #print('analysis_composition_of_labled_sequence_final')
            print('step',13)
            P,C,Sh,U = analysis_composition_of_labled_sequence(plus_sequence,comp_sequence,both_sequecne)

            plus_portions.append(P)
            comp_portions.append(C)
            duble_portions.append(Sh)
            unmask_portions.append(U)
            chr_number_list.append(Chr_id)
            organism_name.append(file_name_spilit[0])
            length_list.append(length)
            p2c_list.append(round(P/C,4))

            #print(f'finish chr{Chr_number,Chr_ID}')
            Chr_sequence_refseq = 0
            plus,comp,both = 0,0,0
            Plus_strand,Complement_strand,Both = 0,0,0
            merged_sequence_plus=0
            plus_sequence = 0
            merged_sequence_comp = 0
            comp_sequence = 0
            merged_sequence_both = 0
            both_sequecne = 0
my_dic = {'organism':organism_name,'chr_number_list':chr_number_list,'length':length_list,
          'plus_portion':plus_portions,'minis_portion':comp_portions,'shered_portion':duble_portions,
          'unmasked_portion':unmask_portions,'P2C':p2c_list}
import pandas as pd
my_datafram = pd.DataFrame(my_dic)
#my_datafram.to_excel(r'D:\PhD\chargaf_genome_test_project\Human\GRCh38.p14\human.xlsx')
my_datafram.to_excel(rf'M:\genome_analysis\100_genome\result\{file_name_spilit[0]}.xlsx')

divition_P_to_C = []
for i in range(len(plus_portions)):
    divition_P_to_C.append(plus_portions[i]/comp_portions[i])
import matplotlib.pyplot as plt2
plt2.plot(divition_P_to_C)
plt2.plot(duble_portions)
plt2.plot(unmask_portions)
plt2.show()
print(divition_P_to_C)
print(plus_portions)
print(comp_portions)


