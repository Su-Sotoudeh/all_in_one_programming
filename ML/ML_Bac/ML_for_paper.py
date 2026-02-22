import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
import joblib

# =====================================================
# 1) LOAD CSV
# =====================================================
csv_path = r"C:\Genomics\results\Genes_chardaff_parameters\Whole_csv_file.csv"
df = pd.read_csv(csv_path, low_memory=False)

# Columns to use
cols_needed = ["A2T_P", "P2C_W", "Length(M)", "A2T_annot"]
for col in cols_needed:
    df[col] = pd.to_numeric(df[col], errors="coerce")

df = df.dropna(subset=cols_needed)
print("Rows after cleaning:", len(df))
if len(df) == 0:
    raise ValueError("Dataframe empty after cleaning; check your CSV.")

# Extract arrays
A2T_P = df["A2T_P"].values
P2C_W = df["P2C_W"].values
Length_M = df["Length(M)"].values
y_real = df["A2T_annot"].values

# =====================================================
# 2) CENTER TARGET: ΔA2T = A2T - 1
# =====================================================
y_centered = y_real - 1

# =====================================================
# 3) CALCULATE FORMULA PREDICTION (centered)
# =====================================================
formula_raw = ((A2T_P * P2C_W) + 1) / (A2T_P + P2C_W)
formula_centered = formula_raw - 1  # center same way

# =====================================================
# 4) SCALE FEATURES
# =====================================================
X = np.column_stack([P2C_W, A2T_P, Length_M])
scaler_X = StandardScaler()
X_scaled = scaler_X.fit_transform(X)

# =====================================================
# 5) TRAIN ML MODEL
# =====================================================
model = RandomForestRegressor(n_estimators=200, random_state=42)
model.fit(X_scaled, y_centered)  # use centered target

# Predict
y_pred_centered = model.predict(X_scaled)

# =====================================================
# 6) METRICS FUNCTION
# =====================================================
def print_metrics(y_true, y_pred, label):
    mse = mean_squared_error(y_true, y_pred)
    rmse = np.sqrt(mse)
    mae = mean_absolute_error(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)
    nrmse_mean = rmse / np.mean(np.abs(y_true))
    nrmse_range = rmse / (np.max(y_true) - np.min(y_true))
    print(f"\n===== {label} =====")
    print(f"R²            : {r2:.4f}")
    print(f"MAE           : {mae:.6f}")
    print(f"RMSE          : {rmse:.6f}")
    print(f"NRMSE (mean)  : {nrmse_mean:.6f}")
    print(f"NRMSE (range) : {nrmse_range:.6f}")

# =====================================================
# 7) PRINT METRICS
# =====================================================
print_metrics(y_centered, formula_centered, "FORMULA (Centered)")
print_metrics(y_centered, y_pred_centered, "ML MODEL (Centered)")

# =====================================================
# 8) PLOTS
# =====================================================
# Scatter
plt.figure(figsize=(7, 7))
plt.scatter(y_centered, formula_centered, alpha=0.4, label="Formula", color="orange")
plt.scatter(y_centered, y_pred_centered, alpha=0.4, label="ML Model", color="green")
low = min(y_centered.min(), formula_centered.min(), y_pred_centered.min())
high = max(y_centered.max(), formula_centered.max(), y_pred_centered.max())
plt.plot([low, high], [low, high], "r--", label="1:1 line")
plt.xlabel("Real ΔA2T")
plt.ylabel("Predicted ΔA2T")
plt.title("Centered Real vs Predictions")
plt.legend()
plt.tight_layout()
plt.show()

# KDE
plt.figure(figsize=(8, 5))
sns.kdeplot(y_centered, label="Real", fill=True)
sns.kdeplot(formula_centered, label="Formula", fill=True)
sns.kdeplot(y_pred_centered, label="ML", fill=True)
plt.title("Distribution Comparison: Real vs Formula vs ML")
plt.xlabel("ΔA2T values")
plt.tight_layout()
plt.show()

# Boxplot
plt.figure(figsize=(6, 5))
sns.boxplot(data=pd.DataFrame({
    "Real": y_centered,
    "Formula": formula_centered,
    "ML": y_pred_centered
}))
plt.title("Boxplot Comparison: Real vs Formula vs ML")
plt.tight_layout()
plt.show()

# =====================================================
# 9) SAVE RESULTS
# =====================================================
output_path = r"C:\Genomics\Chargaff\comparison_A2T_centered.csv"
pd.DataFrame({
    "Real": y_centered,
    "Formula": formula_centered,
    "ML": y_pred_centered
}).to_csv(output_path, index=False)

# Save scaler and model
joblib.dump(scaler_X, r"C:\Genomics\Chargaff\X_scaler.pkl")
joblib.dump(model, r"C:\Genomics\Chargaff\ML_model_A2T_centered.pkl")

print("\nSaved comparison CSV and model successfully.")
