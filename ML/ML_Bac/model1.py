z = [1]*10
print(z)
import pandas as pd
#chordata = pd.read_excel(r'D:\PhD\chargaf_genome_test_project\chordata\refrence_CDS_2025-01-13_06-29-06\results_2025-01-13_06-29-06.xlsx')
#chordata = pd.read_excel(r'D:\Projects\all_in_one_programming\some_result\results_2025-02-05_09-27-33.xlsx')
Bacteria = pd.read_excel(r'D:\PhD\chargaf_genome_test_project\bacteria\CDS_2024-12-29_00-14-30\results_2024-12-29_00-14-30.xlsx')
#Bacteria = Bacteria[(Bacteria["P2C"]<1.05) & (Bacteria["P2C"]>0.95)]

print(Bacteria.head().to_string())
chordata = Bacteria
#chordata = chordata[chordata["Non%"]<50]
#print(chordata.head().to_string())


Y1 = chordata["A2T_P"]
X1 = chordata["P2C"]
Y3 = chordata["A2T_C"]
Y2 = chordata["N_A"]/chordata["N_T"]
X2 = [1]* len(Y2)
M_W_real = chordata["A2T_anotPC_list"]
non_portion = chordata["Non%"]
my_M_W = []
for q,w,e,r,t,m,k in zip(X1,Y1,X2,Y2,non_portion,M_W_real,Y3):
    M_coding = ((q * w) + 1) / (q + w)
    X = q
    Y = w
    A  =( (X/(X+1))  *  (Y/(Y+1))        )   +       (    (1-(X/(X+1))) * ((k/(k+1)))    )
    T  =( (X/(X+1))  *  ((k/(k+1)))  )   +       (    (1-(X/(X+1))) *  (Y/(Y+1))         )
    M = A/T

    M_non_coding = ((e * r) + 1) / (e + r)/100

    M_W = ((100-int(t))*M_coding) + (t*M_non_coding)
    #my_M_W.append(round(M_coding*1000,2))
    my_M_W.append(round(M,3))
    #print(((100-int(t))*M_coding),(t*M_non_coding))
    print(round(M,3),'real',round(m,3),w)



import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from scipy.stats import pearsonr, spearmanr

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
#real_clean = [1]*1000
#pred_clean = real[mask]

print("Original length:", len(real), "Cleaned length:", len(real_clean))



# Errors
'''mse = mean_squared_error(real_clean, pred_clean)        # Mean Squared Error
rmse = np.sqrt(mse)                         # Root MSE
mae = mean_absolute_error(real_clean, pred_clean)       # Mean Absolute Error

# Goodness of fit
r2 = r2_score(real_clean, pred_clean)                   # R² score
pearson_corr, p_value = pearsonr(real_clean, pred_clean) # Pearson correlation
spearman_corr, _ = spearmanr(real_clean, pred_clean)    # Spearman correlation

print("MSE:", mse)
print("RMSE:", rmse)
print("MAE:", mae)
print("R²:", r2)
print("Pearson r:", pearson_corr, " (p=", p_value,")")
print("Spearman ρ:", spearman_corr)


import matplotlib.pyplot as plt

plt.scatter(real_clean, pred_clean, alpha=0.7)
plt.plot([min(real_clean), max(real_clean)], [min(real_clean), max(real_clean)], 'r--')  # 1:1 line
plt.xlabel("Real data")
plt.ylabel("Predicted data")
plt.title("Model Validation")
plt.show()'''




#new method withot R2
import numpy as np
from sklearn.metrics import mean_squared_error, mean_absolute_error
from scipy.stats import pearsonr, spearmanr

def concordance_correlation_coefficient(y_true, y_pred):
    # Lin's CCC
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    mean_true = y_true.mean()
    mean_pred = y_pred.mean()
    var_true = y_true.var(ddof=1)
    var_pred = y_pred.var(ddof=1)
    cov = np.cov(y_true, y_pred, ddof=1)[0,1]
    numerator = 2 * cov
    denominator = var_true + var_pred + (mean_true - mean_pred)**2
    if denominator == 0:
        return np.nan
    return numerator / denominator

def summary_metrics(y_true, y_pred, eps=1e-12):
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    # remove NaN / inf pairs safely
    mask = (~np.isnan(y_true)) & (~np.isnan(y_pred)) & np.isfinite(y_true) & np.isfinite(y_pred)
    y_true = y_true[mask]
    y_pred = y_pred[mask]
    n = len(y_true)
    if n == 0:
        raise ValueError("No valid data points after removing NaN/Inf.")
    mse = mean_squared_error(y_true, y_pred)
    rmse = np.sqrt(mse)
    mae = mean_absolute_error(y_true, y_pred)

    stdev_true = np.std(y_true, ddof=1)
    stdev_pred = np.std(y_pred, ddof=1)

    # safe R^2: compute manually, but detect tiny variance
    ss_res = np.sum((y_true - y_pred)**2)
    ss_tot = np.sum((y_true - y_true.mean())**2)
    if ss_tot <= eps:
        r2 = np.nan  # undefined / meaningless
    else:
        r2 = 1 - ss_res / ss_tot

    # normalized RMSE (by mean and by range)
    nrmse_mean = rmse / (np.abs(y_true.mean()) + eps)
    y_range = y_true.max() - y_true.min()
    nrmse_range = rmse / (y_range + eps)

    # MAPE (avoid division by zero)
    nonzero_mask = np.abs(y_true) > eps
    if nonzero_mask.sum() > 0:
        mape = np.mean(np.abs((y_true[nonzero_mask] - y_pred[nonzero_mask]) / y_true[nonzero_mask])) * 100.0
    else:
        mape = np.nan

    # correlations (Pearson/Spearman) but only if variance > 0
    pearson_r = np.nan
    pearson_p = np.nan
    spearman_rho = np.nan
    try:
        if stdev_true > eps and stdev_pred > eps:
            pearson_r, pearson_p = pearsonr(y_true, y_pred)
            spearman_rho, _ = spearmanr(y_true, y_pred)
    except Exception:
        pass

    ccc = concordance_correlation_coefficient(y_true, y_pred)

    return {
        "n": n,
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "r2": r2,
        "nrmse_mean": nrmse_mean,
        "nrmse_range": nrmse_range,
        "mape_%": mape,
        "pearson_r": pearson_r,
        "pearson_p": pearson_p,
        "spearman_rho": spearman_rho,
        "ccc": ccc,
        "stdev_true": stdev_true,
        "y_range": y_range
    }

# Example for your constant case:
y_true = real_clean        # all ones
y_pred = pred_clean

print(summary_metrics(y_true, y_pred))



import matplotlib.pyplot as plt
def bland_altman_plot(y_true, y_pred):
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    mean_vals = (y_true + y_pred) / 2.0
    diffs = y_true - y_pred
    md = np.mean(diffs)
    sd = np.std(diffs, ddof=1)
    plt.scatter(mean_vals, diffs, alpha=0.5)
    plt.axhline(md, color='red', linestyle='--', label=f'mean diff={md:.4g}')
    plt.axhline(md + 1.96*sd, color='grey', linestyle=':', label='±1.96 SD')
    plt.axhline(md - 1.96*sd, color='grey', linestyle=':')
    plt.xlabel('Mean of true & pred')
    plt.ylabel('True - Pred')
    plt.title('Bland-Altman')
    plt.legend()
    plt.show()

bland_altman_plot(y_true,y_pred)


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Suppose you already have:


# Scatter plot with 1:1 line
plt.figure(figsize=(6,6))
sns.scatterplot(x=y_true, y=y_pred, alpha=0.5)
max_val = max(max(y_true), max(y_pred))
min_val = min(min(y_true), min(y_pred))
plt.plot([min_val, max_val], [min_val, max_val], 'r--', label="1:1 line")

plt.xlabel("True values")
plt.ylabel("Predicted values")
plt.title("Model Predictions vs. True values")
plt.legend()
plt.show()


# Residual plot
plt.figure(figsize=(6,4))
residuals = np.array(y_true) - np.array(y_pred)
sns.scatterplot(x=y_true, y=residuals, alpha=0.5)
plt.axhline(0, color='r', linestyle='--')
plt.xlabel("True values")
plt.ylabel("Residuals (True - Pred)")
plt.title("Residual Plot")
plt.show()


import matplotlib.pyplot as plt
import seaborn as sns

# Assuming y_true and y_pred are numpy arrays or lists

plt.figure(figsize=(8,5))

sns.kdeplot(y_true, label="True values", fill=True, alpha=0.5)
sns.kdeplot(y_pred, label="Predicted values", fill=True, alpha=0.5)

plt.xlabel("Value")
plt.ylabel("Density")
plt.title("Distribution of True vs Predicted Values")
plt.legend()
plt.show()


import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Put into a DataFrame for easier plotting
df = pd.DataFrame({
    "True values": y_true,
    "Predicted values": y_pred
})

plt.figure(figsize=(6,5))
sns.boxplot(data=df, palette="Set2")

plt.title("Boxplot of True vs Predicted Values")
plt.ylabel("Value")
plt.show()

