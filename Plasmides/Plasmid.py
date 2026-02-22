import pandas as pd
from Bio import SeqIO
'''file = open(r"C:\\Users\\rezas\Downloads\plasmid.1.1.genomic.fna\plasmid.1.1.genomic.fna")

seqs = SeqIO.parse(file, "fasta")
for seq in seqs :
    #print(seq.id)
    print(len(seq.seq))'''



import gzip
from Bio import SeqIO

gbff = "C:\Genomics\Genomes\Plasmid\plasmid.4.genomic.gbff.gz"
gff_out = "C:\Genomics\Genomes\Plasmid\plasmid.4.genomic.gff"
df = pd.read_csv(gff_out)


print(df.head())

'''A
with gzip.open(gbff, "rt") as handle, open(gff_out, "w") as out:
    out.write("##gff-version 3\n")

    for record in SeqIO.parse(handle, "genbank"):
        seqid = record.id

        for feature in record.features:
            if feature.type not in ["gene", "CDS", "rRNA", "tRNA"]:
                continue

            start = int(feature.location.start) + 1   # GFF is 1-based
            end = int(feature.location.end)
            strand = "+" if feature.location.strand == 1 else "-"
            phase = "0" if feature.type == "CDS" else "."

            attrs = []

            for key in ["gene", "locus_tag", "product", "protein_id"]:
                if key in feature.qualifiers:
                    attrs.append(f"{key}={feature.qualifiers[key][0]}")

            attributes = ";".join(attrs) if attrs else "."

            out.write(
                f"{seqid}\tGenBank\t{feature.type}\t"
                f"{start}\t{end}\t.\t{strand}\t{phase}\t{attributes}\n")'''

