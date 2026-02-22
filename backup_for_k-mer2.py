

# listing file
import os
import gzip
import pandas as pd
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

save_file = '40_40.csv'
if  os.path.exists(save_file):
    save_file = '40_40.csv'
    new_df = pd.read_csv(save_file,index_col=0)
    plt.figure(figsize=(10, 6))
    ax = sns.heatmap(new_df, annot=False, fmt='2f', cmap='viridis',yticklabels=new_df.index)  # fmt='d' for integer annotations
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,fontsize=10, ha='right') ###
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0,fontsize=8, ha='right')

    plt.title("Heatmap K-mers in 40 organism")
    plt.show()

    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    from matplotlib.patches import ConnectionPatch

    # Set publication-style aesthetics
    sns.set_style('white')  # Clean background
    sns.set_context('paper', font_scale=1.2)  # Adjust font sizes for publication

    # Load the CSV file
    save_file = '40_40.csv'
    new_df = pd.read_csv(save_file, index_col=0)

    # Plot the heatmap
    fig, ax = plt.subplots(figsize=(10, 8))  # Adjusted for taller plot to fit lines/labels
    heatmap = sns.heatmap(new_df, annot=False, fmt='.2f', cmap='viridis', yticklabels=new_df.index,
                          cbar_kws={'label': 'K-mer Frequency'})  # Add colorbar label

    # Simplified x-axis labels: +M, +RC, -M, -RC repeated for each group
    simplified_labels = ['+M', '+RC', '-M', '-RC'] * 3
    ax.set_xticklabels(simplified_labels, rotation=45, fontsize=10, ha='right')

    # Style y-axis labels
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8, ha='right')

    # Axis labels
    ax.set_xlabel('K-mer Variations', fontsize=12)
    ax.set_ylabel('Organisms', fontsize=12)

    # Title (for draft; remove if using figure caption in paper)
    ax.set_title("Heatmap of K-mers in 40 Organisms", fontsize=14)

    # Add vertical lines to separate the 3 groups (after columns 4 and 8)
    ax.axvline(x=4, color='white', linewidth=2)
    ax.axvline(x=8, color='white', linewidth=2)

    # Define group positions
    group_starts = [0.5, 5, 9]  # Start of each group
    group_ends = [3.5, 7, 11]   # End of each group
    group_centers = [1.5, 5.5, 9.5]  # Center of each group
    group_names = ["AAATG", "AAAA", "GGGA"]

    # Add dendrogram-like lines to group the four columns
    for i, (start, end, center) in enumerate(zip(group_starts, group_ends, group_centers)):
        # Independent horizontal line spanning the four columns (one per group), with small gaps to avoid looking connected
        ax.hlines(y=len(new_df) + 2, xmin=start + 0.05, xmax=end - 0.05, color='black', linewidth=1.5)

        # Short vertical line up from center of horizontal to K-mer name
        ax.vlines(x=center, ymin=len(new_df) + 2, ymax=len(new_df) + 2.5, color='black', linewidth=1.5)

        # K-mer name above the center
        ax.text(center, len(new_df) + 2.5, group_names[i], ha='center', va='bottom', fontsize=10, fontweight='bold')

        # Add 4 short vertical lines down from the horizontal line to each column center
        for col in range(start, end):
            ax.vlines(x=col + 0.5, ymin=len(new_df) + 1.5, ymax=len(new_df) + 2, color='black', linewidth=1.5)

    # Add curved arrows to show relations: +M <-> -RC and +RC <-> -M for each group
    for offset in [0, 4, 8]:
        # +M (col offset+0) <-> -RC (col offset+3)
        con1 = ConnectionPatch(xyA=(offset + 0.5, len(new_df) + 1), xyB=(offset + 3.5, len(new_df) + 1),
                               coordsA='data', coordsB='data', axesA=ax, axesB=ax,
                               arrowstyle='<->,head_length=0.5,head_width=0.5', color='black', lw=1.5,
                               connectionstyle='arc3,rad=0.3')
        ax.add_patch(con1)

        # +RC (col offset+1) <-> -M (col offset+2)
        con2 = ConnectionPatch(xyA=(offset + 1.5, len(new_df) + 1), xyB=(offset + 2.5, len(new_df) + 1),
                               coordsA='data', coordsB='data', axesA=ax, axesB=ax,
                               arrowstyle='<->,head_length=0.5,head_width=0.5', color='black', lw=1.5,
                               connectionstyle='arc3,rad=0.3')
        ax.add_patch(con2)

    # Adjust y-limits to make room for lines and labels
    ax.set_ylim(0, len(new_df) + 3)

    # Use tight_layout for automatic centering and spacing
    plt.tight_layout(pad=2.0)  # Add padding to center the figure

    # Save for publication
    plt.savefig('heatmap_publication.pdf', dpi=300, bbox_inches='tight', format='pdf')  # Vector format
    plt.savefig('heatmap_publication.png', dpi=300, bbox_inches='tight')  # High-res raster backup

    plt.show()





else:
    one_Trabite = '/media/su/EC3A11743A113D50/output/' # pc 5118
    #Bioinf = '/media/su/Bioinf/output/Test/' # pc 5118
    #path  = '/media/hpc-sotoudeh/ADATA_HD720/output/Human/' #hpc
    #path = r'D:\PhD\genomes\Mus_19\\' #laptop
    path = one_Trabite
    organisms = os.listdir(path)



    organism_counter = 0

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
        path += organism +'/'
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

        print(fna)
        print(out)
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
        df.to_csv('40_40.csv')
        #print('procces finished')


        # Create the heatmap
        plt.figure(figsize=(10, 6))
        sns.heatmap(df, annot=True, fmt='2f', cmap='viridis',yticklabels=df.index)  # fmt='d' for integer annotations
        plt.title("Heatmap of Your DataFrame")
        plt.show()
