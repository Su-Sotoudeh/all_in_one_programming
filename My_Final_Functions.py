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
        print('opening_fasta_file finished')
        records.append([str(seq_record.id),str(l),str(seq_record.seq)])
    return records