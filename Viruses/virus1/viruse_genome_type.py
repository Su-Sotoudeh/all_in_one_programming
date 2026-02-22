from Bio import Entrez

# Provide your email to NCBI (this is required)
Entrez.email = "reza.sotoudeh.van@gmail.com"

def get_genome_type(genome_id):
    try:
        # Fetch the genome data from NCBI
        handle = Entrez.efetch(db="nucleotide", id=genome_id, rettype="gb", retmode="text")
        record = handle.read()
        handle.close()

        # Extract relevant information
        if "dsDNA" in record:
            return "dsDNA"
        elif "ssDNA" in record:
            return "ssDNA"
        elif "dsRNA" in record:
            return "dsRNA"
        elif "ssRNA" in record:
            return "ssRNA"
        else:
            return "Unknown"
    except Exception as e:
        return f"Error: {str(e)}"

# List of genome IDs
genome_ids = ["HG425345.1", "DQ150593.1"]

# Fetch genome types for all IDs
genome_types = {genome_id: get_genome_type(genome_id) for genome_id in genome_ids}

# Print the results
for genome_id, genome_type in genome_types.items():
    print(f"Genome ID: {genome_id}, Genome Type: {genome_type}")
