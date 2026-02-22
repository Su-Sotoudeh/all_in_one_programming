import pandas as pd
import os
import glob


z = open("bad_files.txt","r")
print(z)
'''organims_list = ["Manis pentadactyla","Vulpes lagopus","Ursus arctos",
                "Vulpes vulpes","Manis javanica","Sarcophilus harrisii",
                 "Notamacropus eugenii","Canis lupus baileyi"
                 ,"Piliocolobus tephrosceles","Piliocolobus tephrosceles","Piliocolobus tephrosceles","Lagenorhynchus albirostris"
                 ,"Phacochoerus africanus","Saccopteryx leptura","Saccopteryx leptura","Balaenoptera ricei","Mesoplodon densirostris","Pseudorca crassidens"
                 ,"Cervus canadensis","Cervus canadensis","Myotis daubentonii","Myotis daubentonii","Diceros bicornis minor","Phyllostomus discolor"

                 ]'''

path = ("/media/su/Bioinf/Genomes_2025_3/Genomes_2025_3_7/chordata_annotated_complete_chr_genome_gff_cds_2025/ncbi_dataset/data")
genome_data_list = os.listdir(path)
TSV_file = pd.read_csv("mammals_ref_anot_2025_8_17.tsv", sep="\t")
organims_list = list(TSV_file["Organism Name"])
def extract_inforamtion_TSV(organism_name,TSV_data_fram):
    row = TSV_data_fram[TSV_data_fram["Organism Name"] == organism_name]
    Organism_Name =row.iloc[0]["Organism Name"].replace(' ',"_")

    Assembly_Accession = row.iloc[0]["Assembly Accession"]

    Taxonomic_ID = row.iloc[0]["Organism Taxonomic ID"]

    #print(Organism_Name,Assembly_Accession,Taxonomic_ID)
    destination_folder = os.path.join('/media/su/EC3A11743A113D50/output', Organism_Name)

    return  Organism_Name,Assembly_Accession,Taxonomic_ID,destination_folder
def check_organism_in_database(Assembly_Accession):
    T_F = True
    if Assembly_Accession in genome_data_list:
        print("we have organism data")
        T_F = True
    else:
        print("we do not have organism data")
        T_F = False
    return T_F
def check_folder_existence(organism_folder_path_in_output):
    A = 0
    if os.path.exists(organism_folder_path_in_output):
        print(f"⚠️ Genome already exists: {organism_folder_path_in_output}")
        print("len files inside organism folder", len(os.listdir(organism_folder_path_in_output)))
        print(os.listdir(organism_folder_path_in_output))
        A = True
    else:

        os.mkdir(organism_folder_path_in_output)
        print(f"📁 Creating folder: {organism_folder_path_in_output}")
        A = False
    return A
def teest_and_copy_gz_file(Organism_name,Assembly_Accession,genome_file_path_in_database,destination_folder) :

    gff_file =  glob.glob(genome_file_path_in_database+"/"+Assembly_Accession+"/*.gff.gz")
    fna_file = glob.glob(genome_file_path_in_database+"/"+Assembly_Accession+"/G*.fna.gz")
    print(genome_file_path_in_database+"/"+Assembly_Accession)
    print(len(gff_file))
    print(len(fna_file))
    zip_genome_file = fna_file[0]
    zip_gff_file = gff_file[0]
    import gzip
    try:
        with gzip.open(zip_genome_file, 'rt') as f:
            # Read the first line (or a few lines to test)
            for _ in range(3):
                print(f.readline().strip())
        print("✅ File is valid and readable.")
        import os
        import shutil
        import gzip



        # === User Inputs ===
        source_file = zip_genome_file       # Path to your .fna.gz file
        destination_root = '/media/su/EC3A11743A113D50/output'              # Where to create the new folder
        new_folder_name = Organism_name                 # Name for the new folder

        # === Step 3: Copy .gz file to the new folder ===
        copied_file_path = os.path.join(destination_folder, os.path.basename(zip_genome_file))
        copied_gff_file_path = os.path.join(destination_folder, os.path.basename(zip_gff_file))
        #basename extrac file name
        shutil.copy2(zip_genome_file, copied_file_path)
        shutil.copy2(zip_gff_file, copied_gff_file_path)
        print(f"📄 Copied: {copied_file_path}")
        print(f"📄 Copied: {copied_gff_file_path}")
        # === Step 4: Extract the .gz file ===
        extracted_file_path = copied_file_path[:-3]  # Remove '.gz'
        extracted_gff_file_path = copied_gff_file_path[:-3]
        try:
            with gzip.open(copied_file_path, 'rb') as f_in:
                with open(extracted_file_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                    os.remove(copied_file_path)
                    print("remove zip file")
            print(f"✅ Extracted to: {extracted_file_path}")
            print(f"✅ Extracted to: {extracted_gff_file_path}")
        except Exception as e:
            print(f"❌ Extraction failed: {e}")
            os.remove(copied_file_path)
            os.remove(destination_root)
            z.write(Organism_name)
            return False


    except Exception as e:
        print("❌ Failed to read the file:", e)


for genome in organims_list:

    organ_Name,Assembly_Accession,Taxonomic_ID,destination_folder = extract_inforamtion_TSV(genome,TSV_file)

    print(organ_Name)
    if check_organism_in_database(Assembly_Accession) == True:


        True_or_False = check_folder_existence(destination_folder)
        if True_or_False == False: #folder creat now and sould continue to analysis
            print("RepeatMasker")
            if len(os.listdir(destination_folder)) == 0:
                A = teest_and_copy_gz_file(organ_Name,Assembly_Accession,path,destination_folder)
                if A == False:
                    continue
            if "G" in os.listdir(destination_folder)[0] and os.listdir(destination_folder)[0].endswith(".fna") :#G is for GCF_ddbfgbfbf
                print("genome_dividing")
                genome_file_name = (os.listdir(destination_folder)[0])
                from RepeatMasker.Genome_file_dividing import Whole_genome_to_chr
                Whole_genome_to_chr(destination_folder+"/"+genome_file_name, destination_folder)
        if True_or_False == True:
            print(destination_folder,print(len(destination_folder)))
            if len(os.listdir(destination_folder)) == 0:
                A = teest_and_copy_gz_file(organ_Name,Assembly_Accession,path,destination_folder)
                if A == False:
                    continue
            if "G" in os.listdir(destination_folder)[0] and os.listdir(destination_folder)[0].endswith(".fna") :#G is for GCF_ddbfgbfbf
                print("genome_dividing")
                genome_file_name = (os.listdir(destination_folder)[0])
                from RepeatMasker.Genome_file_dividing import Whole_genome_to_chr
                Whole_genome_to_chr(destination_folder+"/"+genome_file_name, destination_folder)

                #os.remove(destination_folder+"/"+genome_file_name)
        if "G" not in os.listdir(destination_folder)[0]:
            from RepeatMasker.Linux_in_python_RepeatMasker import RepetMasker
            try:
                RepetMasker(destination_folder,str(Taxonomic_ID))
            except Exception as e:
                print(e)
        if "G" not in os.listdir(destination_folder)[0]:
            print("annotation")
            from RepeatMasker.eukaryotes_for_chargaff import all_element_annotation_in_plus_and_minis
            all_element_annotation_in_plus_and_minis(destination_folder,organ_Name)
    else:
        print("Please download its data")