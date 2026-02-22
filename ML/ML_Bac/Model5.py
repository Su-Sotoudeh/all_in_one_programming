import random

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from scipy.stats import pearsonr, spearmanr


#Bacteria = pd.read_excel(r'C:\Genomics\Chargaff\coding\all_in_one_programming\ML\ML_Bac\CDS_2024-12-29_00-14-30\results_2024-12-29_00-14-30.xlsx')


#Bacteria = Bacteria[(Bacteria["P2C"]<1.05) & (Bacteria["P2C"]>0.95)]




#chordata = chordata[chordata["Non%"]<50]
#print(chordata.head().to_string())
#Bacteria = pd.read_csv(r"C:\Genomics\results\Genes_chardaff_parameters\Bacterial_2025-09-23_20-16-24.csv")
#chordata = pd.read_csv(r"C:\Genomics\results\Genes_chardaff_parameters\Chordata.csv")
#Viride = pd.read_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Viride_on_saveing_time_resuts.csv')
#insect = pd.read_csv(r"C:\Genomics\results\Genes_chardaff_parameters\Arthropoda_500m.csv")
#Fungi = pd.read_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Fungisresults_2025-09-17_00-53-54.csv')
#Archia = pd.read_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Archaea_2025-11-16_20-58-46.csv')
#Chloro_plast = pd.read_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Chlroplast_2025-11-18_22-41-13.csv')
#Mito = pd.read_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Mito_2025-11-18_19-16-56.csv')
#df = pd.concat([chordata,Viride,insect,Fungi,Archia,Bacteria,Chloro_plast,Mito],axis=0,ignore_index=True)

#df =df[(df['Chr_ID'].str.startswith('NC'))&(df['Wh_ano']>10)]


'''maxraj_portion =   df['Int_ano'] +  df['UTR_ano']  + df['CDS_ano'] + df['share_ano']


Y = (df["Int_P_A2T"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_P_A2T"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_P_A2T"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_P_A2T"]*(df['share_ano']/maxraj_portion))
X = (df["Int_P2C"]   * (df['Int_ano']/maxraj_portion)) + (df["UTR_P2C"]  *(df['UTR_ano']/maxraj_portion)) +(df["CDS_P2C"]  *(df['CDS_ano']/maxraj_portion)) + (df["Share_P2C"]  *(df['share_ano']/maxraj_portion))




L = df["Length(M)"]
#  (df['CDS_P_A2T'] * (df['CDS_P%']/(df['CDS_P%']+df['CDS_C%']))) + (df['CDS_C_A2T'] * (df['CDS_C%']/(df['CDS_P%']+df['CDS_C%'])))
M_W_real = (df["Int_A2T_ano"] * (df['Int_ano']/maxraj_portion)) + (df["UTR_A2T_ano"]*(df['UTR_ano']/maxraj_portion)) +(df["CDS_A2T_ano"]*(df['CDS_ano']/maxraj_portion)) + (df["Share_A2T_ano"]*(df['share_ano']/maxraj_portion))

df['A2T_P'] = Y
df['P2C_W'] = X #annotated
df['A2T_annot'] = M_W_real
#df.to_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Whole_data_fram.csv',index=False)'''
df = pd.read_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Whole_csv_file.csv')
my_M_W = []
real_A2T_W = []
predict_A2T_W = []
for q,w,m,l in zip(df['P2C_W'],df['A2T_P'],df['A2T_annot'],df["Length(M)"]):
    try:
        A2T_p = ((q * w) + 1) / (q + w)
    except ZeroDivisionError:
        continue

    print('gene_P2C',q,'gene_P_A2T',w,'gene_A2T_ano',m,'Length(M)',l)
    A_W_predict =   A2T_p/(A2T_p+1) * l
    T_W_predict =  (1- A2T_p/(A2T_p+1)) * l
    A2T_W_predected = A_W_predict/T_W_predict

    A_W_real = A2T_W_predected/(A2T_W_predected+1) *l
    T_W_real = (1 - A_W_real) * l
    A2T_W_real = A_W_real/T_W_real

    #print('L',l,'A_W_real',(m/(m+1)*l),'A_W_predict',A_W_predict)
    #print('A2T_W_real',m,'A_T_W_predict',A_W_predict/T_W_predict)


    #A_portion = m/(m+1)
    #T_portion = 1 - A_portion
    #Extended_A = (A_portion-T_portion)
    real_A2T_W.append(m)


    #A_portionp = A2T_W_predected/(A2T_W_predected+1)
    #T_portionp = 1 - A_portionp
    #Extended_Ap = (A_portionp-T_portionp)
    predict_A2T_W.append(A2T_p)



# Example data
real = real_A2T_W
pred = predict_A2T_W

# Convert to numpy arrays
real = np.array(real, dtype=float)
pred = np.array(pred, dtype=float)

# Remove NaN (and Inf) pairs
mask = ~np.isnan(real) & ~np.isnan(pred) & ~np.isinf(real) & ~np.isinf(pred)
real_clean = real[mask]
pred_clean = pred[mask]




# -----------------------
# Example: your real & predicted data
# Replace with your own arrays/lists
# -----------------------
y_true = np.array(real)  # replace 'real' with your list/array
y_pred = np.array(pred)  # replace 'pred' with your list/array

# Remove NaNs
mask = ~np.isnan(y_true) & ~np.isnan(y_pred)
y_true, y_pred = y_true[mask], y_pred[mask]

# -----------------------
# Metrics
# -----------------------
mse = mean_squared_error(y_true, y_pred)
rmse = np.sqrt(mse)
mae = mean_absolute_error(y_true, y_pred)
r2 = r2_score(y_true, y_pred)

pearson_r, pearson_p = pearsonr(y_true, y_pred)
spearman_r, spearman_p = spearmanr(y_true, y_pred)

# CCC
def concordance_correlation_coefficient(y_true, y_pred):
    mean_true, mean_pred = np.mean(y_true), np.mean(y_pred)
    var_true, var_pred = np.var(y_true), np.var(y_pred)
    cov = np.mean((y_true - mean_true) * (y_pred - mean_pred))
    return (2 * cov) / (var_true + var_pred + (mean_true - mean_pred) ** 2)

ccc = concordance_correlation_coefficient(y_true, y_pred)

# Normalized errors
nrmse_mean = rmse / np.mean(y_true)
nrmse_range = rmse / (np.max(y_true) - np.min(y_true))
mape = np.mean(np.abs((y_true - y_pred) / y_true)) * 100



import numpy as np

def print_regression_metrics_simple(results):
    """
    Simple formatted output without colors - uses actual results values
    """
    print("\n" + "="*60)
    print("REGRESSION MODEL EVALUATION".center(60))
    print("="*60)

    print(f"\nDataset: {results['n']:,} samples")
    print(f"True values - Std: {results['Std(True)']:.4f}, Range: {results['Range(True)']:.4f}")

    print(f"\n{'METRIC':<20} {'VALUE':<15} {'INTERPRETATION'}")
    print("-"*60)

    # Use actual values from the results dictionary
    metrics_info = [
        ("R²", results["R²"], "↑ 0.7+ excellent, 0.5+ good, <0.3 poor"),
        ("Pearson r", results["Pearson r"], "↑ 0.7+ strong, 0.5+ moderate, <0.3 weak"),
        ("Spearman ρ", results["Spearman ρ"], "↑ 0.7+ strong, 0.5+ moderate, <0.3 weak"),
        ("CCC", results["CCC"], "↑ 0.7+ excellent, 0.5+ good, <0.3 poor"),
        ("MSE", results["MSE"], "↓ Lower is better"),
        ("RMSE", results["RMSE"], "↓ Lower is better"),
        ("MAE", results["MAE"], "↓ Lower is better"),
        ("MAPE %", results["MAPE %"], "↓ <1% excellent, <5% good, >10% poor"),
        ("NRMSE_mean", results["NRMSE_mean"], "↓ Lower is better"),
        ("NRMSE_range", results["NRMSE_range"], "↓ Lower is better"),
    ]

    for metric, value, interpretation in metrics_info:
        print(f"{metric:<20} {value:<15.6f} {interpretation}")

# Your results data
# Summary
results = {
    "n": len(y_true),
    "MSE": mse,
    "RMSE": rmse,
    "MAE": mae,
    "R²": r2,
    "Pearson r": pearson_r,
    "Pearson p": pearson_p,
    "Spearman ρ": spearman_r,
    "Spearman p": spearman_p,
    "CCC": ccc,
    "NRMSE_mean": nrmse_mean,
    "NRMSE_range": nrmse_range,
    "MAPE %": mape,
    "Std(True)": np.std(y_true),
    "Range(True)": np.ptp(y_true),
}

# Usage - now it will use your actual metric values
print_regression_metrics_simple(results)


# -----------------------
# Plots
# -----------------------

# 1. Scatter plot
plt.figure(figsize=(6,6))
plt.scatter(y_true, y_pred, alpha=0.5)
plt.plot([y_true.min(), y_true.max()], [y_true.min(), y_true.max()], 'r--', label="1:1 line")
plt.xlabel("True values")
plt.ylabel("Predicted values")
plt.title("Model Predictions vs True values")
plt.legend()
plt.show()

# 2. Residual plot
residuals = y_pred - y_true
plt.figure(figsize=(6,4))
plt.scatter(y_true, residuals, alpha=0.5)
plt.axhline(0, color='red', linestyle='--')
plt.xlabel("True values")
plt.ylabel("Residuals")
plt.title("Residual Plot")
plt.show()

# 3. KDE distribution
plt.figure(figsize=(8,5))
sns.kdeplot(y_true, label="True values", fill=True, alpha=0.5)
sns.kdeplot(y_pred, label="Predicted values", fill=True, alpha=0.5)
plt.xlabel("Value")
plt.ylabel("Density")
plt.title("Distribution of True vs Predicted Values")
plt.legend()
plt.show()

# 4. Boxplot
df = pd.DataFrame({"True values": y_true, "Predicted values": y_pred})
plt.figure(figsize=(6,5))
sns.boxplot(data=df, palette="Set2")
plt.title("Boxplot of True vs Predicted Values")
plt.ylabel("Value")
plt.show()
