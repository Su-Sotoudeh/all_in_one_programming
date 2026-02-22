import joblib

model = joblib.load('my_first_GC_model.joblib')

result = model.predict([[0.40,0.60,1000]])
print(result)
import pandas as pd

dic = pd.DataFrame([{'GC':30,
                    'AT':70,
                    'Length':1000}])


print(model.predict(dic))