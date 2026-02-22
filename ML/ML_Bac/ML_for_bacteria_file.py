import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_selection import mutual_info_regression
from sklearn.metrics import r2_score

# Load data
Bacteria = pd.read_csv(r"C:\Genomics\results\Genes_chardaff_parameters\Bacteria_2025-11-11_17-47-44.csv")
Mito = pd.read_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Mito_2025-11-18_19-16-56.csv')
Archia = pd.read_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Archaea_2025-11-16_20-58-46.csv')
Fungi = pd.read_csv(r'C:\Genomics\results\Genes_chardaff_parameters\Fungisresults_2025-09-17_00-53-54.csv')
insect = pd.read_csv(r"C:\Genomics\results\Genes_chardaff_parameters\Arthropoda_500m.csv")
chordata = pd.read_csv(r"C:\Genomics\results\Genes_chardaff_parameters\Chordata.csv")
Viride_plant = pd.read_csv(r"C:\Genomics\results\Genes_chardaff_parameters\Viride_on_saveing_time_resuts.csv")
df = Viride_plant
#df = df['Wh_ano']
#df = df[((df['CDS_P2C']>0.99)&(df['CDS_P2C']<1.01))&(df['CDS_P_A2T']>1.2)]
#print(df['CDS_P_A2T'],(1/(df['CDS_C_A2T'])))

#df['Semilarity_A2T_PandC'] = (df['CDS_P_A2T']*df['CDS_P2C'])/((1/(df['CDS_C_A2T'])*(1/df['CDS_P2C']))) #* df['CDS_P2C']
#df['Semilarity_A2T_PandC'] = ((df['CDS_P_A2T'] * df['CDS_P2C']) + 1) / (df['CDS_P_A2T'] + df['CDS_P2C'])
df['Semilarity_A2T_PandC'] = ((df['CDS_P_A2T'] * df['CDS_P%'])  / ((1/df['CDS_C_A2T'])*df['CDS_C%']))

#for i,j in zip(df['CDS_P_A2T'],(1/(df['CDS_C_A2T']))):
#    print(i/j)
import matplotlib.pyplot as plt
plt.scatter(df['Semilarity_A2T_PandC']/(df['Semilarity_A2T_PandC'].mean()), df['CDS_A2T_ano'],s=0.5,color='red')
#plt.scatter(df['CDS_P2C'], df['A2T_W'],s=0.5,color='red')
plt.xlabel('Semilarity_A2T_PandC')
plt.ylabel('CDS_A2T_ano')
plt.xlim(0,2)
plt.show()
print(df['Semilarity_A2T_PandC'])

'''q= (df['CDS_P_A2T'] * (df['CDS_P%']/(df['CDS_P%']+df['CDS_C%']))) + (df['CDS_C_A2T'] * (df['CDS_C%']/(df['CDS_P%']+df['CDS_C%'])))
w= df['CDS_P2C']
A2T_p = ((q * w) + 1) / (q + w)
A_W_predict =   A2T_p/(A2T_p+1)
T_W_predict =  (1- A2T_p/(A2T_p+1))
A2T_W_predected = A_W_predict/T_W_predict
df['M'] = A2T_W_predected * df["Length(M)"]

# Remove non-numeric columns
df = df.select_dtypes(include=["float64", "int64"])

# Split
X = df.drop(columns=["A2T_W"])
y = df["A2T_W"] *df["Length(M)"]

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Random Forest
model = RandomForestRegressor(n_estimators=500, random_state=42)
model.fit(X_train, y_train)

# Feature importance
importances = pd.Series(model.feature_importances_, index=X.columns)
importances_sorted = importances.sort_values(ascending=False)

#print(importances_sorted.head(20))

# Mutual information
mi = mutual_info_regression(X, y)
mi_series = pd.Series(mi, index=X.columns).sort_values(ascending=False)
#print(mi_series.head(20))'''
