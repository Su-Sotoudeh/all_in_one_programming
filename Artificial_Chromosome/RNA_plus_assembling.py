import matplotlib.pyplot as plt
import pandas as pd
import random
import numpy as np

df = pd.read_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\Viruses\Viruses_Refrence_Datafram.csv')
'''df = df[df['Molecule_type'] == 'ssRNA(+)']
original_df = df.copy()
#df = df[df['CDS_P_A2T']>1.15]
plt.hist(df['CDS_P_A2T'])
plt.show()
print(df.shape)
print(df.head().to_string())


results = []



# Initialize dictionary to store all results
results_dict = {
    'run': [],
    'A2T_range': [],
    'P2C_range': [],
    'filtered_shape': [],
    'selected_rows': [],
    'remaining_rows': [],
    'A2T_P': [],
    'G2C_P': [],
    'A2T_C': [],
    'G2C_C': [],
    'A2T_anno': [],
    'G2C_anno': [],
    'length_A_sel': [],
    'length_T_sel': [],
    'length_C_sel': [],
    'length_G_sel': [],
    'length_A_rem': [],
    'length_T_rem': [],
    'length_C_rem': [],
    'length_G_rem': [],
    'selected_total': [],
    'remaining_total': [],
    'All_annotation_length': [],
    'A2T_P-A2T_C': []
}


import random
for run in range(6000):
    print(f"\n=== RUN {run+1} ===")


    # Load your distribution data ONCE before the loop
    #a = pd.read_csv(r"C:\Genomics\results\Genes_chardaff_parameters\Bacteria_2025-11-11_17-47-44.csv")
    #A2T_distribution = a['CDS_P_A2T'].dropna().values  # Get values, remove NaNs

    # Then inside your loop:
    #A2T_range = np.random.choice(A2T_distribution)  # Random sample from your distribution
    #print(f"Selection threshold: {A2T_range:.3f}")




    A2T_range = random.randint(60,160)
    A2T_range = A2T_range/100

    # Load your distribution data ONCE before the loop
    #a = pd.read_csv(r"C:\Genomics\results\Genes_chardaff_parameters\Bacteria_2025-11-11_17-47-44.csv")
    #P2C_distribution = a['CDS_P2C'].dropna().values  # Get values, remove NaNs

    # Then inside your loop:
    #randome_P2C = np.random.choice(P2C_distribution)  # Random sample from your distribution
    randome_P2C = random.randint(45,55) /100

    print(f"Selection threshold: {randome_P2C:.3f}")

    df_filtered = original_df[(original_df['CDS_P_A2T']>A2T_range)&(original_df['CDS_P_A2T']<(A2T_range+0.4))].copy()
    #print(df_filtered['CDS_P_A2T'])
    #print(df_filtered.shape)

    # Skip if no rows in this range
    if df_filtered.empty:
        print("⚠️ No rows in this range, skipping...")

        # Append NaN values to dictionary
        results_dict['run'].append(run + 1)
        results_dict['A2T_range'].append(A2T_range)
        results_dict['P2C_range'].append(randome_P2C)
        results_dict['filtered_shape'].append((0, 0))
        results_dict['selected_rows'].append(0)
        results_dict['remaining_rows'].append(0)
        results_dict['A2T_P'].append(float('nan'))
        results_dict['G2C_P'].append(float('nan'))
        results_dict['A2T_C'].append(float('nan'))
        results_dict['G2C_C'].append(float('nan'))
        results_dict['A2T_anno'].append(float('nan'))
        results_dict['G2C_anno'].append(float('nan'))
        results_dict['length_A_sel'].append(float('nan'))
        results_dict['length_T_sel'].append(float('nan'))
        results_dict['length_C_sel'].append(float('nan'))
        results_dict['length_G_sel'].append(float('nan'))
        results_dict['length_A_rem'].append(float('nan'))
        results_dict['length_T_rem'].append(float('nan'))
        results_dict['length_C_rem'].append(float('nan'))
        results_dict['length_G_rem'].append(float('nan'))
        results_dict['selected_total'].append(float('nan'))
        results_dict['remaining_total'].append(float('nan'))
        results_dict['All_annotation_length'].append(float('nan'))
        results_dict['A2T_P-A2T_C'].append(float('nan'))
        continue

    annotation_poerion = df_filtered['CDS_P%']/100

    L = df_filtered['Length(M)']

    A = df_filtered['CDS_P_A_ano']/100
    T = df_filtered['CDS_P_T_ano']/100
    C = df_filtered['CDS_P_C_ano']/100
    G = df_filtered['CDS_P_G_ano']/100

    An = L * annotation_poerion * A
    Tn = L * annotation_poerion * T
    Cn = L * annotation_poerion * C
    Gn = L * annotation_poerion * G
    Annotation_sequeence_plus = [sum(An), sum(Tn), sum(Cn), sum(Gn)]
    All_annotation_length_val = sum(Annotation_sequeence_plus)




    # Start fresh with filtered data for each run
    remaining_df = df_filtered.copy()
    selected_rows_list = []

    # Initialize counters for SELECTED group
    length_A_sel = length_T_sel = length_C_sel = length_G_sel = 0

    # Continue until we exceed the threshold OR no rows left



    # Then in your while loop:
    while len(remaining_df) > 0 and (length_A_sel + length_T_sel + length_C_sel + length_G_sel) <= All_annotation_length_val * randome_P2C:

        # Randomly select ONE row by Name/ID
        selected_row = remaining_df.sample(n=1)
        selected_id = selected_row['Name'].iloc[0]  # Assuming 'Name' is your ID column

        # Get the row data
        row = selected_row.iloc[0]

        # Add to selected counters
        length_A_sel += row['Length(M)'] * round(row['CDS_P%']/100, 3) * round(row['CDS_P_A_ano']/100, 3)
        length_T_sel += row['Length(M)'] * round(row['CDS_P%']/100, 3) * round(row['CDS_P_T_ano']/100, 3)
        length_C_sel += row['Length(M)'] * round(row['CDS_P%']/100, 3) * round(row['CDS_P_C_ano']/100, 3)
        length_G_sel += row['Length(M)'] * round(row['CDS_P%']/100, 3) * round(row['CDS_P_G_ano']/100, 3)

        # Store selected row info
        selected_rows_list.append(selected_row)

        # REMOVE this row from remaining dataframe
        remaining_df = remaining_df[remaining_df['Name'] != selected_id]

    #print(f"Selected ID: {selected_id}, Running total: {length_A_sel+length_T_sel+length_C_sel+length_G_sel:.2f}")

    # Create DataFrames
    selected_df = pd.concat(selected_rows_list, ignore_index=True) if selected_rows_list else pd.DataFrame()

    # Calculate totals for REMAINING group
    if len(remaining_df) > 0:
        length_A_rem = sum(remaining_df['Length(M)'] * round(remaining_df['CDS_P%']/100, 3) * round(remaining_df['CDS_P_A_ano']/100, 3))
        length_T_rem = sum(remaining_df['Length(M)'] * round(remaining_df['CDS_P%']/100, 3) * round(remaining_df['CDS_P_T_ano']/100, 3))
        length_C_rem = sum(remaining_df['Length(M)'] * round(remaining_df['CDS_P%']/100, 3) * round(remaining_df['CDS_P_C_ano']/100, 3))
        length_G_rem = sum(remaining_df['Length(M)'] * round(remaining_df['CDS_P%']/100, 3) * round(remaining_df['CDS_P_G_ano']/100, 3))
    else:
        length_A_rem = length_T_rem = length_C_rem = length_G_rem = 0

    selected_total = length_A_sel + length_T_sel + length_C_sel + length_G_sel
    remaining_total = length_A_rem + length_T_rem + length_C_rem + length_G_rem

    # Safe division calculations
    if length_T_sel > 0:
        a2t_p_val = round(length_A_sel/length_T_sel, 3)
    else:
        a2t_p_val = float('nan')

    if length_C_sel > 0:
        g2c_p_val = round(length_G_sel/length_C_sel, 3)
    else:
        g2c_p_val = float('nan')

    if length_A_rem > 0:
        a2t_c_val = round(length_T_rem/length_A_rem, 3)
    else:
        a2t_c_val = float('nan')

    if length_G_rem > 0:
        g2c_c_val = round(length_C_rem/length_G_rem, 3)
    else:
        g2c_c_val = float('nan')

    if (length_T_sel + length_A_rem) > 0:
        a2t_anno_val = round((length_A_sel+length_T_rem)/(length_T_sel+length_A_rem), 3)
    else:
        a2t_anno_val = float('nan')

    if (length_C_sel + length_G_rem) > 0:
        g2c_anno_val = round((length_G_sel+length_C_rem)/(length_C_sel+length_G_rem), 3)
    else:
        g2c_anno_val = float('nan')

    A2T_P_A2T_C = a2t_p_val - a2t_c_val

    # Append ALL results to dictionary
    results_dict['run'].append(run + 1)
    results_dict['A2T_range'].append(A2T_range)
    results_dict['P2C_range'].append(randome_P2C)
    results_dict['filtered_shape'].append(df_filtered.shape)
    results_dict['selected_rows'].append(len(selected_df))
    results_dict['remaining_rows'].append(len(remaining_df))
    results_dict['A2T_P'].append(a2t_p_val)
    results_dict['G2C_P'].append(g2c_p_val)
    results_dict['A2T_C'].append(a2t_c_val)
    results_dict['G2C_C'].append(g2c_c_val)
    results_dict['A2T_anno'].append(a2t_anno_val)
    results_dict['G2C_anno'].append(g2c_anno_val)
    results_dict['length_A_sel'].append(length_A_sel)
    results_dict['length_T_sel'].append(length_T_sel)
    results_dict['length_C_sel'].append(length_C_sel)
    results_dict['length_G_sel'].append(length_G_sel)
    results_dict['length_A_rem'].append(length_A_rem)
    results_dict['length_T_rem'].append(length_T_rem)
    results_dict['length_C_rem'].append(length_C_rem)
    results_dict['length_G_rem'].append(length_G_rem)
    results_dict['selected_total'].append(selected_total)
    results_dict['remaining_total'].append(remaining_total)
    results_dict['All_annotation_length'].append(All_annotation_length_val)
    results_dict['A2T_P-A2T_C'].append(A2T_P_A2T_C)

    #print(f"Selected: {len(selected_df)} rows, A/T: {a2t_p_val if not pd.isna(a2t_p_val) else 'NaN'}")
    #print(f"Remaining: {len(remaining_df)} rows, T/A: {a2t_c_val if not pd.isna(a2t_c_val) else 'NaN'}")

# Convert dictionary to DataFrame
results_df = pd.DataFrame(results_dict)
results_df = results_df.sort_values( ['A2T_P-A2T_C'])
print("\n" + "="*80)
print("RESULTS DATAFRAME")
print("="*80)
print(results_df.head(20).to_string())
print(f"\nTotal runs recorded: {len(results_df)}")
print(results_df.describe().to_string())

# Save to CSV
results_df.to_csv('simulation_results.csv', index=False)
print("\nResults saved to 'simulation_results.csv'")'''
results_df = pd.read_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\Artificial_Chromosome\simulation_results.csv')
# Plot using DataFrame columns
x = np.arange(len(results_df))
plt.figure(figsize=(12, 6))
plt.scatter(x, results_df['A2T_P'], color='r', s=1, label='A/T +')
plt.scatter(x, results_df['A2T_C'], color='g', s=1, label='A/T -')
plt.scatter(x, results_df['A2T_anno'], color='black', s=1, label='A/T Anno')
plt.legend()
plt.title('A/T Ratios Across 6000 Runs')
plt.xlabel('Run')
plt.ylabel('Ratio')
plt.show()

# Additional plots
plt.figure(figsize=(12, 6))
plt.scatter(x, results_df['G2C_P'], color='r', s=1, label='G/C +')
plt.scatter(x, results_df['G2C_C'], color='g', s=1, label='C/G -')
plt.scatter(x, results_df['G2C_anno'], color='black', s=1, label='G/C Anno')
plt.legend()
plt.title('G/C Ratios Across 6000 Runs')
plt.xlabel('Run')
plt.ylabel('Ratio')
plt.show()

# Histogram of A2T ratios
'''plt.figure(figsize=(10, 6))
plt.hist(results_df['A2T_P'].dropna(), bins=30, alpha=0.5, label='A/T +', color='r')
plt.hist(results_df['A2T_C'].dropna(), bins=30, alpha=0.5, label='A/T -', color='g')
plt.hist(results_df['A2T_anno'].dropna(), bins=30, alpha=0.5, label='A/T Anno', color='black')
plt.legend()
plt.title('Distribution of A/T Ratios')
plt.xlabel('Ratio')
plt.ylabel('Frequency')
plt.show()'''