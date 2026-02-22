import gzip
import pandas as pd
from Bio import SeqIO

def opening_fasta_file(seq_path):
    fasta_path = seq_path
    if fasta_path.endswith('gz'):
        fasta_file = gzip.open(fasta_path,'rt')
    else:
        fasta_file = fasta_path
    records = []
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(seq_record.seq)
        l  = len(sequence)/1000000
        print('sequence id :',seq_record.id,' len chr :',l,'Mb')
        #print('opening_fasta_file finished')
        records.append(str(seq_record.seq))
    return records




def Plus_minus_lableing(sequence, gff_path):

    sequence = str(sequence)
    #print(len(sequence))
    sequence = sequence.lower()  # ATCG all
    sequence = sequence.replace('a','1')
    sequence = sequence.replace('t','2')
    sequence = sequence.replace('c','3')
    sequence = sequence.replace('g','4')

    sequence_plus = sequence
    sequence_plus = list(sequence_plus) #ATCG all



    #lableing genes***************************************************************************
    if gff_path.endswith('.gz'):
        gff = gzip.open(gff_path)
    else:
        gff = open(gff_path)
    df = pd.read_csv(gff,sep='\t',comment='#')
    df.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    #print(df.head().to_string())
    counter2 = 0
    filtered_df = df[ df['type'].str.startswith('gene')]
    for annotation in filtered_df.to_numpy():

        counter2 += 1
        start = int(annotation[3])
        end = int(annotation[4])
        if end - start == 0:
            print('len was 0')
            continue
        if annotation[6] == '+':
            gene_sequence = sequence_plus[(start - 1):end]
            gene_sequence = ''.join(gene_sequence)
            gene_sequence = gene_sequence.replace('1','A')
            gene_sequence = gene_sequence.replace('2','T')
            gene_sequence = gene_sequence.replace('3','C')
            gene_sequence = gene_sequence.replace('4','G')
            gene_sequence = gene_sequence.replace('a','Q')
            gene_sequence = gene_sequence.replace('t','W')
            gene_sequence = gene_sequence.replace('c','E')
            gene_sequence = gene_sequence.replace('g','R')
            sequence_plus[(start - 1):end] = gene_sequence

        if annotation[6] == '-':
            gene_sequence = sequence_plus[(start - 1):end]
            gene_sequence = ''.join(gene_sequence)
            gene_sequence = gene_sequence.replace('1','a')
            gene_sequence = gene_sequence.replace('2','t')
            gene_sequence = gene_sequence.replace('3','c')
            gene_sequence = gene_sequence.replace('4','g')
            gene_sequence = gene_sequence.replace('A','Q')
            gene_sequence = gene_sequence.replace('T','W')
            gene_sequence = gene_sequence.replace('C','E')
            gene_sequence = gene_sequence.replace('G','R')
            sequence_plus[(start - 1):end] = gene_sequence

    gff.close()

    #print('gene pass',len(sequence_plus))
    sequence_plus = ''.join(sequence_plus)
    return sequence_plus







genome_path = r"C:\Genomics\Genomes\Aureus_test3\ncbi_dataset\data\GCF_000013425.1\GCF_000013425.1_ASM1342v1_genomic.fna"
gff_path = r"C:\Genomics\Genomes\Aureus_test3\ncbi_dataset\data\GCF_000013425.1\genomic.gff"


sequence = opening_fasta_file(genome_path)[0]
#print(sequence[:1000])
annntated_sequence = Plus_minus_lableing(sequence, gff_path)
#print(annntated_sequence[10000:15000])

'''group2 = ['A','T','C','G','a','t','c','g','Q','W','E','R']
for i in group2:
    print(annntated_sequence.count(i))'''
