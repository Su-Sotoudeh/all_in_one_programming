import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report

# Load CSV
df = pd.read_csv("genes_gc.csv")

# Features and labels
X = df[["GC", "AT", "Length"]]
y = df["label"]

# Train/test split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Model
model = RandomForestClassifier(n_estimators=200)

# Train
model.fit(X_train, y_train)

# Predict
pred = model.predict(X_test)

# Accuracy
print("Accuracy:", accuracy_score(y_test, pred))
print(classification_report(y_test, pred))

import matplotlib.pyplot as plt

plt.scatter(df["GC"], df["AT"], c=(df["label"]=="high"))
plt.xlabel("GC")
plt.ylabel("AT")
plt.title("High vs Low GC Genes")
plt.show()

import joblib

joblib.dump(model, 'my_first_GC_model.joblib')