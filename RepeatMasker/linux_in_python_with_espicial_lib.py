import os
import subprocess

# Directory containing the FASTA files
input_dir = os.path.abspath("/media/su/EC3A11743A113D50/output_sumer2025/Arabidopsis")

# Get all .fna files in the input directory
fasta_files = [f for f in os.listdir(input_dir) if f.endswith('.fna')]

# Process each FASTA file
for fasta_file in fasta_files:
    input_file = os.path.join(input_dir, fasta_file)
    base_name, _ = os.path.splitext(fasta_file)
    out_file = os.path.join(input_dir, f"{base_name}.fna.out")


    # Skip if .out file already exists
    if os.path.exists(out_file):
        print(f"Skipping {fasta_file} — .out file already exists.")
        continue

    print(f"Running RepeatMasker on {fasta_file}...")

    command = [
        "/usr/local/RepeatMasker/RepeatMasker",

        "-no_is", "-nolow", "-norna",
        "-pa", "7", "-engine", "rmblast", '-lib', "/media/su/Bioinf/Dfam38/RepeatMasker.lib",
        input_file
    ]

    subprocess.run(command, check=True)
for fasta_file in fasta_files:
    input_file = os.path.join(input_dir, fasta_file)
    base_name, _ = os.path.splitext(fasta_file)
    # Remove .masked file if it exists
    masked_file = os.path.join(input_dir, f"{base_name}.fna.masked")

    if os.path.exists(masked_file):
        os.remove(masked_file)
        print(f"Removed: {masked_file}")

    # Remove .cat.gz file if it exists
    cat_gz_file = os.path.join(input_dir, f"{base_name}.fna.cat.gz")
    if os.path.exists(cat_gz_file):
        os.remove(cat_gz_file)
        print(f"Removed: {cat_gz_file}")
