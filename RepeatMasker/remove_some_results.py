#********************************************************csv remove!!!!!! ****************************************

#input_dir = os.listdir("/media/hpc-sotoudeh/ADATA_HD720/output/")
path = '/media/su/EC3A11743A113D50/output/'
input_dir = os.listdir(path)

print(input_dir)
for organism in input_dir:
    print(organism)
    fasta_files = [f for f in os.listdir(f"/media/su/EC3A11743A113D50/output/{organism}") if f.endswith('.cat.gz') or f.endswith('.masked') or f.endswith('.csv')]
    if len(fasta_files) == 0:
        continue

    print(fasta_files)

    for file in fasta_files:
        file_path = f'/media/su/EC3A11743A113D50/output/{organism}/{file}'
        print(f'deleted{file}')
        os.remove(file_path)



