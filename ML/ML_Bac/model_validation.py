
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from scipy.stats import pearsonr, spearmanr


Bacteria = pd.read_excel(r'D:\PhD\chargaf_genome_test_project\bacteria\CDS_2024-12-29_00-14-30\results_2024-12-29_00-14-30.xlsx')
#Bacteria = Bacteria[(Bacteria["P2C"]<1.05) & (Bacteria["P2C"]>0.95)]

print(Bacteria.head().to_string())
chordata = Bacteria
#chordata = chordata[chordata["Non%"]<50]
#print(chordata.head().to_string())
chordata = pd.read_excel(r'D:\PhD\chargaf_genome_test_project\chordata\refrence_CDS_2025-01-13_06-29-06\results_2025-01-13_06-29-06.xlsx')
insect = pd.read_excel(r"D:\Projects\all_in_one_programming\some_result\insect\refrence_CDS_2025-02-02_05-57-17\results_2025-02-02_05-57-17.xlsx")
chordata = insect

Y1 = chordata["A2T_P"]
X1 = chordata["P2C"]
Y3 = chordata["A2T_C"]
Y2 = chordata["N_A"]/chordata["N_T"]

X2 = [1]* len(Y2)
M_W_real = chordata["A2T_W"].round(2)
non_portion = chordata["Non%"]
my_M_W = []
for q,w,e,r,t,m,k in zip(X1,Y1,X2,Y2,non_portion,M_W_real,Y3):
    M_coding = ((q * w) + 1) / (q + w)
    #X = q
    #Y = w
    #A  =( (X/(X+1))  *  (Y/(Y+1))        )   +       (    (1-(X/(X+1))) * ((k/(k+1)))    )
    #T  =( (X/(X+1))  *  ((k/(k+1)))  )   +       (    (1-(X/(X+1))) *  (Y/(Y+1))         )
    #M = A/T

    M_non_coding = ((e * r) + 1) / (e + r)
    M_W = ((100-t)*M_coding) + (t*M_non_coding)
    M_W =   round( M_W/100,2)
    print(round(100-t,2),round(M_coding,2) , round(t,2),M_non_coding,"A2T_W:",round(m,2),"my_A2T_W:",M_W)

    my_M_W.append(M_W)
    #print(((100-int(t))*M_coding),(t*M_non_coding))
    #print("my",round(M_W,3),'real',round(m,3))



# Example data
real = M_W_real
pred = my_M_W



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
