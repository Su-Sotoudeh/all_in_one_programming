# -----------------------------
# Full RandomForest ML for A2T (3 features)
# -----------------------------

import warnings
warnings.filterwarnings("ignore")

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from scipy.stats import pearsonr, spearmanr
import joblib

# -----------------------------
# 1) LOAD DATA
# -----------------------------
'''file_paths = {
    "Bacteria": r"C:\Genomics\results\Genes_chardaff_parameters\Bacteria_2025-11-11_17-47-44.csv",
    "chordata": r"C:\Genomics\results\Genes_chardaff_parameters\Chordata.csv",
    "Viride": r'C:\Genomics\results\Genes_chardaff_parameters\Viride_on_saveing_time_resuts.csv',
    "insect": r"C:\Genomics\results\Genes_chardaff_parameters\Arthropoda_500m.csv",
    "Fungi": r'C:\Genomics\results\Genes_chardaff_parameters\Fungisresults_2025-09-17_00-53-54.csv',
    "Archia": r'C:\Genomics\results\Genes_chardaff_parameters\Archaea_2025-11-16_20-58-46.csv',
    "Chloro_plast": r'C:\Genomics\results\Genes_chardaff_parameters\Chlroplast_2025-11-18_22-41-13.csv',
    "Mito": r'C:\Genomics\results\Genes_chardaff_parameters\Mito_2025-11-18_19-16-56.csv'
}

dfs = []
for name, path in file_paths.items():
    try:
        d = pd.read_csv(path)
        dfs.append(d)
    except Exception:
        print(f"Warning: couldn't load {name} from {path} (skipped).")

if len(dfs) == 0:
    raise RuntimeError("No input dataframes found. Check file paths.")

df = pd.concat(dfs, axis=0, ignore_index=True)
print(f"Total rows after concat: {len(df):,}")

# -----------------------------
# 2) CALCULATE WEIGHTED FEATURES
# -----------------------------
maxraj_portion = df['Int_ano'] + df['UTR_ano'] + df['CDS_ano'] + df['share_ano']
mask_valid_portion = maxraj_portion > 0
df = df[mask_valid_portion].copy()
maxraj_portion = maxraj_portion[mask_valid_portion]

Y = (df["Int_P_A2T"] * (df['Int_ano'] / maxraj_portion)) + \
    (df["UTR_P_A2T"] * (df['UTR_ano'] / maxraj_portion)) + \
    (df["CDS_P_A2T"] * (df['CDS_ano'] / maxraj_portion)) + \
    (df["Share_P_A2T"] * (df['share_ano'] / maxraj_portion))

X = (df["Int_P2C"] * (df['Int_ano'] / maxraj_portion)) + \
    (df["UTR_P2C"] * (df['UTR_ano'] / maxraj_portion)) + \
    (df["CDS_P2C"] * (df['CDS_ano'] / maxraj_portion)) + \
    (df["Share_P2C"] * (df['share_ano'] / maxraj_portion))

L = df["Length(M)"]

M_W_real = (df["Int_A2T_ano"] * (df['Int_ano'] / maxraj_portion)) + \
           (df["UTR_A2T_ano"] * (df['UTR_ano'] / maxraj_portion)) + \
           (df["CDS_A2T_ano"] * (df['CDS_ano'] / maxraj_portion)) + \
           (df["Share_A2T_ano"] * (df['share_ano'] / maxraj_portion))'''

# -----------------------------
# 3) BUILD ML DATAFRAME
# -----------------------------

df = pd.read_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Whole_csv_file.csv')

ml_df = pd.DataFrame({
    "P2C_w": df['P2C_W'],
    "P_A2T_w": df['A2T_P'],
    "Length_M": df["Length(M)"],
    "Real_A2T_ano": df['A2T_annot']
})
ml_df = ml_df.replace([np.inf, -np.inf], np.nan).dropna()
ml_df = ml_df.reset_index(drop=True)
print(f"ML dataset rows (after dropna): {len(ml_df):,}")

X_ml = ml_df[["P2C_w", "P_A2T_w", "Length_M"]].values
y_ml = ml_df["Real_A2T_ano"].values.reshape(-1, 1)

# Scale target
y_scaler = StandardScaler()
y_ml_scaled = y_scaler.fit_transform(y_ml).ravel()

# Save scaler
joblib.dump(y_scaler, r"C:\Genomics\Chargaff\y_scaler_A2T.pkl")
print("Target scaler saved!")

# -----------------------------
# 4) TRAIN-TEST SPLIT
# -----------------------------
X_train, X_test, y_train_scaled, y_test_scaled = train_test_split(
    X_ml, y_ml_scaled, test_size=0.2, random_state=42
)
y_test_original = y_scaler.inverse_transform(y_test_scaled.reshape(-1,1)).ravel()

# -----------------------------
# 5) RANDOM FOREST + GRIDSEARCH
# -----------------------------
param_grid = {
    "n_estimators": [200, 400],
    "max_depth": [None, 10, 20],
    "min_samples_split": [2, 5],
}

rf = RandomForestRegressor(random_state=42, n_jobs=-1)
grid = GridSearchCV(rf, param_grid, cv=3, scoring="r2", n_jobs=-1, verbose=0)
grid.fit(X_train, y_train_scaled)
best_rf = grid.best_estimator_
print("GridSearch best params:", grid.best_params_)

# Save model
joblib.dump(best_rf, r"C:\Genomics\Chargaff\ML_model_A2T.pkl")
print("Model saved successfully!")

# -----------------------------
# 6) PREDICT & INVERSE-SCALE
# -----------------------------
y_pred_scaled = best_rf.predict(X_test)
y_pred_original = y_scaler.inverse_transform(y_pred_scaled.reshape(-1,1)).ravel()

# -----------------------------
# 7) METRICS
# -----------------------------
def concordance_correlation_coefficient(y_true, y_pred):
    mean_true, mean_pred = np.mean(y_true), np.mean(y_pred)
    var_true, var_pred = np.var(y_true), np.var(y_pred)
    cov = np.mean((y_true - mean_true) * (y_pred - mean_pred))
    return (2 * cov) / (var_true + var_pred + (mean_true - mean_pred) ** 2)

y_true = y_test_original
y_pred = y_pred_original

mse = mean_squared_error(y_true, y_pred)
rmse = np.sqrt(mse)
mae = mean_absolute_error(y_true, y_pred)
r2 = r2_score(y_true, y_pred)
pearson_r, _ = pearsonr(y_true, y_pred)
spearman_r, _ = spearmanr(y_true, y_pred)
ccc = concordance_correlation_coefficient(y_true, y_pred)

print(f"R2: {r2:.4f}, RMSE: {rmse:.4f}, MAE: {mae:.4f}, Pearson r: {pearson_r:.4f}, CCC: {ccc:.4f}")

# -----------------------------
# 8) NEW SAMPLE PREDICTION
# -----------------------------
new_samples = pd.DataFrame({
    "P2C_w": [0.98, 1.02, 1.00],
    "P_A2T_w": [0.8, 1.0, 1.2],
    "Length_M": [1000, 1500, 2000]  # realistic values
})

pred_scaled = best_rf.predict(new_samples.values)
pred_original = y_scaler.inverse_transform(pred_scaled.reshape(-1,1)).ravel()
new_samples["Predicted_A2T"] = pred_original
print("\nPredictions for new samples:")
print(new_samples)
