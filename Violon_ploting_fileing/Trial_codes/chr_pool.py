import time
from Bio import  SeqIO
import os

def prepar_chromosomes(all_genomes_folder_path,organ_group):
    print('inalysis_started')
    start_time = time.time()
    current_time = datetime.now()
    start_run = current_time.strftime("%Y-%m-%d_%H-%M-%S")
    print(start_run)
    organism_list = os.listdir(all_genomes_folder_path)
    organism_number = len(organism_list)
    counter = 0
    import pandas as pd
    saved_file = pd.read_csv('/media/su/EC3A11743A113D50/genome_analysis/fungi/on_saveing_time_resuts.csv')
    the_last_organ = saved_file.tail(8)['Name'].tolist()
    print('the_last_organ',the_last_organ)
    saved_chr_ids = saved_file['Chr_ID'].tolist()
    saved_names = saved_file['Name'].tolist()
    for i in the_last_organ:
        saved_names.remove(i)
    for organism in organism_list:
        counter += 1
        print(organism_number,counter)
        if organism in saved_names  :
            print(organism,'before analysed')
            continue
        if organism.endswith(".json") or organism.endswith(".jsonl"):
            continue
        files = os.listdir(all_genomes_folder_path + '/' + organism)
        fna = [i for i in files if i.endswith('.fna.gz') and not i.startswith('cds')]
        gff = [k for k in files if k.endswith('.gff.gz')]
        fna.sort()
        gff.sort()
        try:
            fasta_file = gzip.open(all_genomes_folder_path + '/' + organism + '/' + fna[0], 'rt')
        except (IOError, EOFError, IndexError):
            return None
        # Load GFF once per organism to save memory (avoid reloading multiple times per chromosome)
        gff_path_main = all_genomes_folder_path + '/' + organism + '/' + gff[0]
        chromosomes_pool = []
        chr_id_pool = []
        organ_name_pool = []
        gff_path_pool = []
        try:
            for seq_record in SeqIO.parse(fasta_file, "fasta"):


                ch_id = seq_record.id
                if ch_id in saved_chr_ids:
                    print(ch_id," this chromosome already analysed")
                    continue
                if ch_id in saved_chr_ids:
                    print('chr_exit_and analysed')
                    continue
                if 'mitochondrion' in seq_record.description:
                    continue
                seq = str(seq_record.seq)
                l = len(seq) / 1000000
                if l < 0.05:
                    continue
                if "chromosome" not in seq_record.description:
                    continue
                if ch_id.startswith('NW'):
                    continue
                chromosomes_pool.append(seq)
                chr_id_pool.append(ch_id)
                organ_name_pool.append(organism)
                gff_path_pool.append(gff_path_main)
                if len(chromosomes_pool) == 16:
                    universal_Pool = 0
                    while universal_Pool == 0:
                        universal_Pool = 0


