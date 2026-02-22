import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import r2_score, mean_absolute_error
from joblib import Parallel, delayed
import multiprocessing

# ==============================
# CONFIGURATION
# ==============================
csv_path = "G:\CDS_UTR_Intr\Chordata.csv"  # <-- change to your file path
target_col = "A2T_W"  # <-- or G2C_W, CDS_P2C, etc.

# ==============================
# LOAD AND PREPARE DATA
# ==============================
df = pd.read_csv(csv_path)

# Remove non-numeric or identifier columns
non_features = ["Name", "Chr_ID"]
feature_cols = [c for c in df.columns if c not in non_features + [target_col]]

# Drop rows with missing values (optional)
df = df.dropna(subset=[target_col] + feature_cols)

X = df[feature_cols]
y = df[target_col]

# ==============================
# TRAIN / TEST SPLIT
# ==============================
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

# ==============================
# RANDOM FOREST MODEL
# ==============================
n_cores = multiprocessing.cpu_count()  # use all available cores
model = RandomForestRegressor(
    n_estimators=400,
    max_depth=None,
    n_jobs=n_cores,      # parallelize across cores
    random_state=42
)

model.fit(X_train, y_train)

# ==============================
# EVALUATION
# ==============================
y_pred = model.predict(X_test)
r2 = r2_score(y_test, y_pred)
mae = mean_absolute_error(y_test, y_pred)

print(f"R² score: {r2:.4f}")
print(f"Mean Absolute Error: {mae:.4f}")

# ==============================
# FEATURE IMPORTANCE
# ==============================
importances = pd.Series(model.feature_importances_, index=feature_cols)
top_features = importances.sort_values(ascending=False).head(20)
print("\nTop 20 most important features:")
print(top_features)

# Optional: Save importances for later analysis
top_features.to_csv("ML_Chordata.csv")

# ==============================
# CROSS-VALIDATION (Parallel)
# ==============================
cv_scores = cross_val_score(model, X, y, cv=5, scoring='r2', n_jobs=n_cores)
print("\nCross-validation R² scores:", cv_scores)
print("Mean CV R²:", np.mean(cv_scores))



import matplotlib.pyplot as plt
top_features.plot(kind='barh')
plt.title("Top 20 Feature Importances")
plt.show()


