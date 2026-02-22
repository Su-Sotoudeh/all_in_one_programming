path = r"C:\Genomics\Chargaff\coding\all_in_one_programming\some_result\results_2025-02-05_09-27-33.xlsx"
#path = "/media/su/EC3A11743A113D50/genome_analysis/chordata/refrence_CDS_2025-01-13_06-29-06/results_2025-01-13_06-29-06.xlsx"
import pandas as pd
dataframe = pd.read_excel(path)
print(dataframe.head().to_string())

import matplotlib.pyplot as plt
sorted_dataframe = dataframe.sort_values(by="A2T_P")
filtered_dataframe0 = sorted_dataframe[(sorted_dataframe["A2T_W"] > 0.99) & (dataframe["A2T_W"] < 1.01)]
filtered_dataframe = filtered_dataframe0[(filtered_dataframe0["Non%"] < 30) ]
import numpy as np

Y = filtered_dataframe["A2T_P"]
Y2 = filtered_dataframe["A2T_C"]
Y3 = filtered_dataframe["A2T_W"]
Y4 = filtered_dataframe["P2C"]
X = np.arange(len(Y))

plt.scatter(X,Y,color="green")
plt.scatter(X,Y2,color="red")
plt.scatter(X,Y3,color="black")
#plt.scatter(X,Y4,color="blue")
plt.show()