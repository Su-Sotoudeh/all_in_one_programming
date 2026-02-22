
import os
import gzip
import pandas as pd
from Bio import SeqIO


def process_one_organism(organism, all_genomes_folder_path):

    files = os.listdir(os.path.join(all_genomes_folder_path, organism))
    fna = [f for f in files if (f.endswith('.fna.gz') or f.endswith('.fna')) and not f.startswith('cds')]
    gff = [f for f in files if (f.endswith('.gff.gz') or f.endswith('.gff'))]

    if not fna or not gff:
        return None

    fna_path = os.path.join(all_genomes_folder_path, organism, fna[0])
    gff_path = os.path.join(all_genomes_folder_path, organism, gff[0])

    fasta_file = gzip.open(fna_path, 'rt') if fna_path.endswith('.gz') else open(fna_path, 'rt')
    gff_file = gzip.open(gff_path, 'rt') if gff_path.endswith('.gz') else open(gff_path, 'r')

    try:
        gff_df = pd.read_csv(gff_file, sep='\t', comment='#',
                             names=['seqid','source','type','start','end','score','strand','phase','attributes'])
    except Exception:
        return None

    results = []
    for seq_record in SeqIO.parse(fasta_file, "fasta"):

        #if 'plasmid' in seq_record.description or 'mitochondrion' in seq_record.description:
        #    continue
        ch_id = seq_record.id
        #if ch_id.startswith('NW') or "chromosome" not in seq_record.description:
        #    continue

        sequence = str(seq_record.seq).upper()
        filtered_df = gff_df[gff_df['seqid'] == ch_id]
        filtered_df = filtered_df[filtered_df['type'] == 'gene']
        genes_sequence_plus, genes_sequence_minis = '', ''
        plus, minis, P_length, M_length = 0, 0, 0, 0

        for _, row in filtered_df.iterrows():
            start, end = int(row["start"]) - 1, int(row["end"])
            gene_sequence = sequence[start:end]
            length = end - start
            if row["strand"] == "+":
                genes_sequence_plus += gene_sequence
                plus += 1
                P_length += length
            elif row["strand"] == "-":
                genes_sequence_minis += gene_sequence
                minis += 1
                M_length += length

        if plus + minis == 0:
            continue

        total_length = len(genes_sequence_plus + genes_sequence_minis)
        data = {
            'P': len(genes_sequence_plus) / total_length * 100 if total_length > 0 else 0,
            'A_plus': genes_sequence_plus.count('A') / P_length * 100 if P_length > 0 else 0,
            'T_plus': genes_sequence_plus.count('T') / P_length * 100 if P_length > 0 else 0,
            'C_plus': genes_sequence_plus.count('C') / P_length * 100 if P_length > 0 else 0,
            'G_plus': genes_sequence_plus.count('G') / P_length * 100 if P_length > 0 else 0,
            'M': len(genes_sequence_minis) / total_length * 100 if total_length > 0 else 0,
            'A_minus': genes_sequence_minis.count('A') / M_length * 100 if M_length > 0 else 0,
            'T_minus': genes_sequence_minis.count('T') / M_length * 100 if M_length > 0 else 0,
            'C_minus': genes_sequence_minis.count('C') / M_length * 100 if M_length > 0 else 0,
            'G_minus': genes_sequence_minis.count('G') / M_length * 100 if M_length > 0 else 0,
            'Chordata': organism
        }
        results.append(data)

    fasta_file.close()
    gff_file.close()
    #print(pd.DataFrame(results).to_string())
    return pd.DataFrame(results)
















genomes_path = r"C:\Genomics\Genomes\unzip_viruse_all_complete\unzip_viruse_all_complete\ncbi_dataset\data\\"
NCBI_information = r"C:\Genomics\Genomes\unzip_viruse_all_complete\unzip_viruse_all_complete\Viruse_file_information\Viruses_information_NCBI_whole.csv"
ICTV_information =  r"C:\Genomics\Genomes\unzip_viruse_all_complete\unzip_viruse_all_complete\Viruse_file_information\viruses_genome_information.csv"
Genome_Type = r"C:\Genomics\Genomes\unzip_viruse_all_complete\unzip_viruse_all_complete\Viruse_file_information\Viruses_Genome_type.csv"
#ICTV_information = pd.read_csv(ICTV_information)
#NCBI_information = pd.read_csv(NCBI_information)
Genome_Type = pd.read_csv(Genome_Type)
#Assembly_Accession = NCBI_information["Assembly Accession"].to_list()
#Organism_Name_NCBI = NCBI_information["Organism Name"].to_list()
#Organism_Name_ICTV = ICTV_information["Virus name(s)"].to_list()'''

my_viruse_data = os.listdir(genomes_path)
filtered_df = Genome_Type[(Genome_Type['Genome'] == 'ssRNA(+)') & (Genome_Type["Assembly Accession"].isin(my_viruse_data))]
print('Genome_number',filtered_df.shape)

viruses_list = filtered_df["Assembly Accession"].to_list()
V_number = len(viruses_list)
print(V_number)
counter = 0
for viruse in viruses_list:
    if viruse.endswith('.json') or viruse.endswith('.jsonl'):
        continue

    result = Genome_Type.loc[Genome_Type["Assembly Accession"] == viruse, "Genome"]

    if isinstance(result, pd.Series) and not result.empty:
        genome_value = result.iloc[0]
    else:
        genome_value = None  # not found

    if genome_value == 'ssRNA(+)':
        files_name = os.listdir(os.path.join(genomes_path[:-1] + viruse))
        counter +=1
        print(counter,viruse, genome_value,len(files_name))

        #fna = os.path.join(genomes_path[:-1] + viruse + '\\' + files_name[0])
        #gff = os.path.join(genomes_path[:-1] + viruse + '\\' + files_name[1])

        #for record in SeqIO.parse(fna,'fasta'):
        #print(counter,V_number,len(record.seq),record.description[:50],str(record.seq)[:50])
        print(process_one_organism(viruse,genomes_path[:-1]).to_string())
