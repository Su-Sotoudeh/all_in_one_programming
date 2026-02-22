import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
import joblib

# =====================================================
# 1) LOAD CSV
# =====================================================
csv_path = r"C:\Genomics\results\Genes_chardaff_parameters\Whole_csv_file.csv"
df = pd.read_csv(csv_path, low_memory=False)


from sklearn.preprocessing import StandardScaler

X_train = df[['P2C_W', 'A2T_P', 'Length(M)']].values
X_scaler = StandardScaler()
X_scaler.fit(X_train)

# Save it for later
joblib.dump(X_scaler, r"C:\Genomics\Chargaff\x_scaler_A2T.pkl")




# Convert necessary columns to numeric
cols_needed = ["A2T_P", "P2C_W", "Length(M)", "A2T_annot"]
for col in cols_needed:
    df[col] = pd.to_numeric(df[col], errors="coerce")

# Remove rows with missing values
df = df.dropna(subset=cols_needed)

print("Rows after cleaning:", len(df))
if len(df) == 0:
    raise ValueError("ERROR: Dataframe empty after cleaning; check your CSV.")

# Extract values
A2T_P = df['A2T_P'].values
P2C_W = df['P2C_W'].values
Length_M = df['Length(M)'].values
y_real = df['A2T_annot'].values


# =====================================================
# 2) FORMULA PREDICTION
# =====================================================
calculated_A2T_by_formula = ((A2T_P * P2C_W) + 1) / (A2T_P + P2C_W)

# Remove infinite or nan values
mask = np.isfinite(calculated_A2T_by_formula)
A2T_P = A2T_P[mask]
P2C_W = P2C_W[mask]
Length_M = Length_M[mask]
y_real = y_real[mask]
calculated_A2T_by_formula = calculated_A2T_by_formula[mask]


# =====================================================
# 3) LOAD MODEL & SCALER
# =====================================================
model_path = r"C:\Genomics\Chargaff\ML_model_A2T.pkl"
scaler_path = r"C:\Genomics\Chargaff\y_scaler_A2T.pkl"

loaded_model = joblib.load(model_path)
y_scaler = joblib.load(scaler_path)

# ML input
X_ml = np.column_stack([P2C_W, A2T_P, Length_M])

# Model prediction
y_pred_scaled = loaded_model.predict(X_ml)
y_pred_ml = y_scaler.inverse_transform(y_pred_scaled.reshape(-1, 1)).ravel()


# =====================================================
# 4) METRIC FUNCTION
# =====================================================asqweaasd
def print_metrics(y_true, y_pred, label):
    mse = mean_squared_error(y_true, y_pred)
    rmse = np.sqrt(mse)
    mae = mean_absolute_error(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)

    nrmse_mean = rmse / np.mean(y_true)
    nrmse_range = rmse / (np.max(y_true) - np.min(y_true))

    print(f"\n===== {label} =====")
    print(f"R²            : {r2:.4f}")
    print(f"MAE           : {mae:.6f}")
    print(f"RMSE          : {rmse:.6f}")
    print(f"NRMSE (mean)  : {nrmse_mean:.6f}")
    print(f"NRMSE (range) : {nrmse_range:.6f}")


# =====================================================
# 5) PRINT PERFORMANCE
# =====================================================
print_metrics(y_real, calculated_A2T_by_formula, "FORMULA Prediction")
print_metrics(y_real, y_pred_ml, "ML MODEL Prediction")


# =====================================================
# 6) PLOTS
# =====================================================

# ------- Scatter Plot --------
plt.figure(figsize=(7, 7))
plt.scatter(y_real, calculated_A2T_by_formula, alpha=0.4, label="Formula", color="orange")
plt.scatter(y_real, y_pred_ml, alpha=0.4, label="ML Model", color="green")

low = min(y_real.min(), calculated_A2T_by_formula.min(), y_pred_ml.min())
high = max(y_real.max(), calculated_A2T_by_formula.max(), y_pred_ml.max())
plt.plot([low, high], [low, high], "r--", label="1:1 line")

plt.xlabel("Real A2T_annot")
plt.ylabel("Predicted A2T")
plt.title("Real vs Predictions (Formula & ML)")
plt.legend()
plt.xlim(0,3.5)
plt.ylim(0,3.5)
plt.tight_layout()
plt.show()


# ------- KDE Distribution --------
plt.figure(figsize=(8, 5))
sns.kdeplot(calculated_A2T_by_formula, label="Formula", fill=True,color="black")
sns.kdeplot(y_real, label="Real", fill=True,color="green")

#sns.kdeplot(y_pred_ml, label="ML", fill=True,color="blue")
plt.title("Distribution Comparison: Real vs Formula")
plt.xlabel("A2T values")
plt.xlim(0.25,1.75)
plt.tight_layout()
plt.legend()

plt.savefig(r'C:\Genomics\Chargaff\furmula_A2T.png',dpi=900)
plt.show()


# ------- KDE Distribution --------
plt.figure(figsize=(8, 5))
#sns.kdeplot(calculated_A2T_by_formula, label="Formula", fill=True,color="black")


sns.kdeplot(y_pred_ml, label="ML", fill=True,color="blue")
sns.kdeplot(y_real, label="Real", fill=True,color="green")
plt.title("Distribution Comparison: Real vs ML")
plt.xlabel("A2T values")
plt.xlim(0.25,1.75)
plt.tight_layout()
plt.legend()
plt.savefig('C:\Genomics\Chargaff\ML_model_A2T.png',dpi=900)
plt.show()




# ------- Boxplot --------
plt.figure(figsize=(6, 5))
sns.boxplot(data=pd.DataFrame({
    "Real": y_real,
    "Formula": calculated_A2T_by_formula,
    "ML": y_pred_ml
}))
plt.title("Boxplot Comparison: Real vs Formula vs ML")
plt.ylim(0.75,1.25)
plt.tight_layout()
plt.show()


# =====================================================
# 7) SAVE RESULTS
# =====================================================
output_path = r"C:\Genomics\Chargaff\comparison_A2T.csv"
pd.DataFrame({
    "Real": y_real,
    "Formula": calculated_A2T_by_formula,
    "ML": y_pred_ml
}).to_csv(output_path, index=False)

print("\nSaved comparison CSV to:", output_path)
print("✓ DONE — All plots generated successfully.")
