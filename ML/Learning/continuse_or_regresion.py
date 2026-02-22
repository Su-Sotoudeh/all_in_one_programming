from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
import pandas as pd

df = pd.read_csv("genes_gc.csv")
X = df[[ "AT", "Length"]]
y = df["GC"]

# Train/test split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
model = RandomForestRegressor(n_estimators=200)
model.fit(X_train,y_train)
pred = model.predict(X_test)


print('MSE:', mean_squared_error(y_test, pred))
print('R2:', r2_score(y_test, pred))

result = model.predict(pd.DataFrame([{'AT':60,'Length':1000}]))
print(result)

