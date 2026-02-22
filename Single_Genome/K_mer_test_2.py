
#new_df = new_df.drop("Drosophila")
#new_df = new_df.drop("C_elegans")
import os

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import ConnectionPatch
from matplotlib.colorbar import Colorbar


# Load the CSV file
save_file = '../112_60.csv'
if os.path.exists(save_file):
    new_df = pd.read_csv(save_file, index_col=0)
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import ConnectionPatch

    # Set publication-style aesthetics
    sns.set_style('white')  # Clean background
    sns.set_context('paper', font_scale=1.2)  # Adjust font sizes for publication

    # Load the CSV file
      # Update to your new file name if different (e.g., '48_40.csv')
    new_df = pd.read_csv(save_file, index_col=0)

    # Define thresholds (adjust these based on your data)
    min_threshold = 5  # Values below this turn white
    max_threshold = 50  # Values above this turn black

    # Create a copy of the data and apply thresholding
    thresholded_df = new_df.copy()
    thresholded_df[(thresholded_df < min_threshold) | (thresholded_df > max_threshold)] = np.nan

    # Plot the heatmap with thresholding
    fig, ax = plt.subplots(figsize=(20, 8))  # Increased width to accommodate 48 columns
    heatmap = sns.heatmap(thresholded_df, annot=False, fmt='.2f', cmap='viridis', yticklabels=new_df.index,
                          cbar_kws={'label': 'K-mer Frequency'}, vmin=min_threshold, vmax=max_threshold)

    # Set tick locations for all 48 columns
    ax.set_xticks(range(60))  # Set ticks for all 48 columns

    # Simplified x-axis labels: +M, +RC, -M, -RC repeated for each group
    simplified_labels = ['+M', '+RC', '-M', '-RC'] * 15  # 48 labels for 12 groups
    ax.set_xticklabels(simplified_labels, rotation=45, fontsize=6, ha='right')  # Reduced fontsize for readability

    # Style y-axis labels
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=4, ha='right')

    # Axis labels
    ax.set_xlabel('K-mer Variations', fontsize=12)
    ax.set_ylabel('Organisms', fontsize=12)

    # Title placement (centered above the plot)
    ax.text(30, len(new_df) + 4.5, "Heatmap of K-mers in 112 Organisms", ha='center', va='bottom', fontsize=14, fontweight='bold')

    # Add vertical lines to separate the 12 groups (after every 4 columns)
    for i in range(4, 60, 4):
        ax.axvline(x=i, color='white', linewidth=2)

    # Define group positions for 12 groups
    group_starts = list(range(0, 60, 4))  # [0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44]
    group_ends = list(range(4, 64, 4))    # [4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48]
    group_centers = [i + 1.5 for i in group_starts]  # [1.5, 5.5, 9.5, 13.5, 17.5, 21.5, 25.5, 29.5, 33.5, 37.5, 41.5, 45.5]
    kmers = ['AT','GG','CA','AAT','CCT','TGA','ATCG','CAAG','TGAC','AAATG','GCGCG','AATTG','CAAGAA','CGCCGA','TTGCGAA']#,'ATATCGC',]

    # Add dendrogram-like lines to group the four columns for all 12 groups
    for i, (start, end, center) in enumerate(zip(group_starts, group_ends, group_centers)):
        # Independent horizontal line spanning the four columns (one per group)
        ax.hlines(y=len(new_df) + 1, xmin=start + 0.472, xmax=end - 0.472, color='0.3', linewidth=1)

        # Short vertical line up from center of horizontal to K-mer name
        #ax.vlines(x=center, ymin=len(new_df) + 2, ymax=len(new_df) + 2.5, color='black', linewidth=1.5)

        # K-mer name above the center
        ax.text(center+0.5, len(new_df) + 1.5, kmers[i], ha='center', va='bottom', fontsize=8, fontweight='bold',color = '0.3')

        # Add 4 short vertical lines down from the horizontal line to each column center
        for col in range(start, end):
            ax.vlines(x=col + 0.5, ymin=len(new_df) , ymax=len(new_df) + 1, color='0.3', linewidth=1)


    # Adjust y-limits to make room for lines and labels
    ax.set_ylim(0, len(new_df) + 3)

    # Use tight_layout for automatic centering and spacing
    plt.tight_layout(pad=2.0)  # Add padding to center the figure

    # Save for publication
    plt.savefig('heatmap_publication.pdf', dpi=600, bbox_inches='tight', format='pdf')  # Vector format
    plt.savefig('heatmap_publication.png', dpi=600, bbox_inches='tight')  # High-res raster backup
    plt.show()


# listing file
import os
import gzip
import pandas as pd
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

save_file = rf'C:\Genomics\Chargaff\coding\all_in_one_programming\{save_file}'
if  os.path.exists(save_file):
    print('yes the 140_40.csv file exists')


    new_df = pd.read_csv(save_file,index_col=0)
    #sorted_df = new_df.sort_values(by=['AAAA_+M'], ascending=False)
    #print(new_df.to_string())
    #new_df.reset_index(drop=False, inplace=True)
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt

    # ----------------------
    # Your existing dataframe with MultiIndex
    # ----------------------
    # Split columns into MultiIndex: (kmer, feature)
    tuples = [tuple(col.split("_")) for col in new_df.columns]
    new_df.columns = pd.MultiIndex.from_tuples(tuples, names=["kmer", "feature"])

    # ----------------------
    # Step 1. Cluster at k-mer level
    # ----------------------
    # Collapse the 4 feature columns into one value per k-mer (e.g., mean)
    df_grouped = new_df.groupby(level="kmer", axis=1).mean()

    # Run clustermap to get the k-mer clustering order
    g = sns.clustermap(
        df_grouped,
        cmap="coolwarm",
        col_cluster=True,
        row_cluster=False,
        figsize=(14, 8)
    )

    # Extract clustered order of kmers
    ordered_kmers = df_grouped.columns[g.dendrogram_col.reordered_ind]

    plt.close()  # close the temporary clustermap figure

    # ----------------------
    # Step 2. Reorder original dataframe by clustered k-mers
    # ----------------------
    features = ["+M", "+RC", "-M", "-RC"]
    new_df = new_df[[(k, f) for k in ordered_kmers for f in features]]

    # ----------------------
    # Step 3. Plot final heatmap with grouped columns
    # ----------------------
    plt.figure(figsize=(14, 8))
    sns.heatmap(
        new_df,
        cmap="coolwarm",
        annot=False,
        yticklabels=new_df.index
    )

    plt.xticks(rotation=45, ha="right", fontsize=9)
    plt.title("K-mer grouped heatmap (clustered by k-mers)")
    plt.tight_layout()

    plt.show()




    plt.figure(figsize=(11, 6))
    #ax = sns.heatmap(sorted_df, annot= False, fmt='.2f', cmap='coolwarm',yticklabels=sorted_df.index)  # fmt='d' for integer annotations
    #ax.set_xticklabels(ax.get_xticklabels(), rotation=45,fontsize=10, ha='right')
    ###ax.set_yticklabels(ax.get_yticklabels(), rotation=0,fontsize=8, ha='right')
    #plt.title("Heatmap of Your DataFrame")
    tuples = [tuple(col.split("_")) for col in new_df.columns ]
    new_df.columns = pd.MultiIndex.from_tuples(tuples,names=["kmer","feature"])
    g = sns.clustermap(new_df,cmap="coolwarm",col_cluster=True,row_cluster=False,annot = False,figsize=(14, 8))
    plt.tight_layout()
    plt.show()

else:
    one_Trabite = r'C:\Genomics\Chargaff\coding\all_in_one_programming\OUT_CSV_files\\'  # pc 5118
    #Bioinf = '/media/su/Bioinf/output/Test/' # pc 5118
    #path  = '/media/hpc-sotoudeh/ADATA_HD720/output/Human/' #hpc
    #path = r'D:\PhD\genomes\Mus_19\\' #laptop
    path = one_Trabite
    organisms = os.listdir(path)



    organism_counter = 0


    kmers = ['GG','CA','AAT','CCT','ATCG','CAAG','AAATG','GGGGG','AAAAAA','CGGGGA','GCGCGCG','ATATCGC']

    kmers = ['AAATG','AAAA','GGGA']

    data_dict = {}
    for k in kmers:
        for feature in ["+M", "+C", "+R", "+RC"]:
            data_dict[f"{k}_{feature}"] = []



    analysed_organisms = []
    for organism in organisms:


        if organism_counter == 41:
            break
        path = one_Trabite
        print(organism,' analysis')
        path += organism +'\\'
        files = os.listdir(path)

        fna = [i for i in files if i.endswith('.fna')]
        out = [j for j in files if j.endswith('.out')]
        gff = [k for k in files if k.endswith('.gff.gz')]
        csv = [p for p in files if p.endswith('.csv')]
        fna.sort()
        out.sort()
        gff.sort()
        from all_functions import  exon_and_intron_masking
        from all_functions import opening_fasta_file
        from all_functions import nucleotide_analysis_for_exon
        from all_functions import TE_masking_after_exon
        from all_functions import nucleotide_analysis_for_exon_TE
        #print('list_of_chromosome for analysis')

        #print('fna',fna)
        if len(out) == 0:
            print('out',out)

        if len(fna) != len(out) or len(out) == 0:
            continue
        for f,o in zip(fna, out):
            #print('chr:',chr_counter,'  FNA',f,'  OUT:',o,'    GFF:',gff[0])
            if f[:11] != o[:11]:
                print('not equal FNA and OUT')
                print('Eroor')
                import sys
                sys.exit()
        else:
            print('all files are match')
        chr_couner = 0
        organism_counter += 1
        print('start anakysis organism', organism_counter)

        if organism_counter == 41:
            break


        for f,o in zip(fna, out):



            chr_couner+=1

            if chr_couner >= 2:
                con = 0
                break
            print('start chr analysed',chr_couner)
            ch_id = f[:-4]
            #if (ch_id + '.csv') in csv:
            #   print(ch_id,'befor proccessd and have csv result')
            #  continue
            #print('start new analysis',f,gff[0],o,ch_id)
            #optional


            seq = opening_fasta_file(path+f)
            #print('A/T',(seq.count('A')+seq.count("a"))/(seq.count('T')+seq.count("t")),'C/G',(seq.count('C')+seq.count("c"))/(seq.count("G")+seq.count("g")),seq.count('N')+seq.count("n"),len(seq))
            exon_mask_seq = exon_and_intron_masking(opening_fasta_file(path+f),path+gff[0],ch_id)
            #print('exon masking pass')
            #nucleotide_analysis_for_exon(exon_mask_seq)
            exon_TE_masked_sequecne = TE_masking_after_exon(exon_mask_seq,path+o,ch_id)
            #print('TE masking pass')
            list_letters = ['ZXVB','SDFH','QWER','*YUI','OPJK','LM@#','1234',
                            'zxvb','sdfh','qwer','+yui','opjk','lm$%','5678',
                            'ATCG','atcg','~~~~']
            df = nucleotide_analysis_for_exon_TE(exon_TE_masked_sequecne)
            #print(df.to_string(formatters={col: '{:,.2f}'.format for col in df.columns}))

            #return_to_ATCG




            #reverse converting to ATCG and atcg
            '''CDS+UTR+Intron+LINE+SINE+LTR+DNA+
            CDS-UTR-Intron-LINE-SINE-LTR-DNA-Shared_Unannotated_unsequenced'''
            k_mers = ['AAATG','AAAA','GGGA']
            list_k_mers = [["***YI","U*YYY","+++yi","u+yyy"],['****','YYYY','++++','yyyy'],['III*','YUUU','iii+','yuuu']]
            # Suppose we have 2 organisms and 3 k-mers

            for index,l in enumerate(k_mers):
                Mp  = exon_TE_masked_sequecne.count(list_k_mers[index][0])
                RCp = exon_TE_masked_sequecne.count(list_k_mers[index][1])
                Mc  = exon_TE_masked_sequecne.count(list_k_mers[index][2])
                RCc = exon_TE_masked_sequecne.count(list_k_mers[index][3])

                pMp  = round(Mp / (Mp+RCp+Mc+RCc) * 100,2 )
                pRCp = round(RCp/ (Mp+RCp+Mc+RCc) * 100,2 )
                pMc  = round(Mc / (Mp+RCp+Mc+RCc) * 100,2 )
                pRCc = round(RCc/ (Mp+RCp+Mc+RCc) * 100,2 )

                data_dict[f"{l}_+M"].append(pMp)
                data_dict[f"{l}_+RC"].append(pRCp)
                data_dict[f"{l}_-M"].append(pMc)
                data_dict[f"{l}_-RC"].append(pRCc)
                print(l,data_dict)
        analysed_organisms.append(organism)
        print(analysed_organisms)
        print(data_dict)
        df = pd.DataFrame(data_dict, index=analysed_organisms)
        print('df is ...')
        print(df)



        df.index = analysed_organisms
        #df.to_csv('K_mer_heatmap_2_40.csv')
        print(df.to_string(formatters={col: '{:,.0f}'.format for col in df.columns}))

        df.to_csv('140_40.csv')

        #print('procces finished')


        # Create the heatmap
        plt.figure(figsize=(10, 6))
        sns.heatmap(df, annot=True, fmt='2f', cmap='viridis',yticklabels=df.index)  # fmt='d' for integer annotations
        plt.title("Heatmap of Your DataFrame")
        plt.show()
