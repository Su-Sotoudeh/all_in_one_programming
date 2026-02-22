input_file = "/media/su/EC3A11743A113D50/output/Mastomys_coucha/GCF_008632895.1_UCSF_Mcou_1_genomic.fna"
output_dir = "/media/su/EC3A11743A113D50/output/Mastomys_coucha"

def Whole_genome_to_chr(file_path, output_path):
    import gzip
    from Bio import SeqIO
    import os

    # Input file (FASTA format)
    input_file = file_path
    # Output directo
    output_dir = output_path
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    #input_file = gzip.open(input_file, 'rt')
    # Read the input file and split sequences
    for record in SeqIO.parse(input_file, "fasta"):
        #print(record.description)
        # Check if 'chromosome' is in the record description (case-insensitive)
        #print(record.id,len(record.seq))
        if record.id[:2] == "NC" :
            print(record.id)
            #print(record.description)
        if 'scaffold' in record.description or 'NW_' in record.description or '>NW' in record.description or 'patch' in record.description:
            continue
        if "chromosome" in record.description.lower():
            chrom_name = record.id  # Get chromosome name
            output_file = os.path.join(output_dir, f"{chrom_name}.fna")  # Create output file path

            # Write the chromosome sequence to a separate file
            with open(output_file, "w") as out_f:
                SeqIO.write(record, out_f, "fasta")

            print(f"Saved: {output_file}")
        else:
            print(f"Skipped scaffold: {record.id}")
    os.remove(file_path)

#Whole_genome_to_chr(input_file, output_dir)