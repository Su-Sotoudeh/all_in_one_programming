import os
import shutil

# Set your paths
input_dir = "/media/su/EC3A11743A113D50/output/"
output_dir = "/media/su/Bioinf/CSV_files_bioinf_2025_9_30"

# Create the CSV_files folder if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Loop over each organism
for organism in os.listdir(input_dir):
    org_path = os.path.join(input_dir, organism)
    if os.path.isdir(org_path):
        # Create a subfolder in CSV_files
        dest_folder = os.path.join(output_dir, organism)
        os.makedirs(dest_folder, exist_ok=True)

        # Copy only .csv files
        for file in os.listdir(org_path):
            if file.endswith(".csv"):
                src_file = os.path.join(org_path, file)
                dst_file = os.path.join(dest_folder, file)
                shutil.copy2(src_file, dst_file)
