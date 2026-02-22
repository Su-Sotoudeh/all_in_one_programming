#need to edit thay not be sort
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.pyplot import legend

file = pd.read_excel(r'C:\Genomics\Chargaff\coding\all_in_one_programming\ML\ML_Bac\results_2024-12-17_17-52-12\results_2024-12-17_17-52-12.xlsx')

file = pd.read_csv(r'C:\Genomics\Chargaff\coding\all_in_one_programming\whole_datafram\Whole_datafram_2026.csv')
#file = file[(file['Group'] != 'Mito') & (file['Group'] != 'Chloro_plast')]
organims_name = 'Chloro_plast'
file = file[file['Group'] == organims_name]
import numpy as np
'''P_sort = list(file['S2B_2'])
P_sort.sort()
x = np.arange(1,len(P_sort)+1)
x2 = [44.5]*len(file['P%'])
x3 = [10]*len(file['P%'])'''





import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

# Example data
#x = file['P_A_usage_list']
#y = file['C_T_usage_list']
#Z = file['P_T_usage_list']

x = file['CDS_P_A_ano']
y = file['CDS_C_T_ano']
Z = file['CDS_P_T_ano']

# Scatter plot
plt.scatter(x, y,alpha=0.5,color = 'green',s = 0.3)
plt.scatter(x, Z,color = 'black',s=1)
line = np.arange(0,max(max(file['CDS_P_A_ano']),max(file['CDS_C_T_ano'])),0.01)
plt.scatter(line,line,color = 'gray',alpha=0.5,s=0.5)
plt.xlabel('Ap usage')
plt.ylabel('Tc usage')
plt.title(f'{organims_name}_Ap to Tc usage')
plt.show()

x = file['P_T_usage_list']
y = file['C_A_usage_list']
Z = file['P_A_usage_list']
# Scatter plot
plt.scatter(x, y,alpha=0.5,color = 'green',s = 0.3)
plt.scatter(x, Z,color = 'black',s=1)
plt.xlabel('Tp usage')
plt.ylabel('Ac usage')
plt.title('Tp to Ac usage')
plt.show()


x = file['P_G_usage_list']
y = file['C_C_usage_list']
Z = file['P_C_usage_list']
# Scatter plot
plt.scatter(x, y,alpha=0.5,color = 'green',s = 0.3)
plt.scatter(x, Z,color = 'black',s=1)
plt.xlabel('Gp usage')
plt.ylabel('Cc usage')
plt.title('Gp to Cc usage')
plt.show()


x = file['P_C_usage_list']
y = file['C_G_usage_list']
Z = file['P_G_usage_list']
# Scatter plot
plt.scatter(x, y,alpha=0.5,color = 'green',s = 0.3)
plt.scatter(x, Z,color = 'black',s=1)
plt.xlabel('Cp usage')
plt.ylabel('Gc usage')
plt.title('Cp to Gc usage')
plt.show()





# Pearson Correlation Coefficient
pearson_corr, _ = pearsonr(x, y)
print(f"Pearson Correlation Coefficient: {pearson_corr:.2f}")

# Spearman Correlation Coefficient
spearman_corr, _ = spearmanr(x, y)
print(f"Spearman Correlation Coefficient: {spearman_corr:.2f}")




from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

# Example data
x = file[['P_A_usage_list']]

y = file['C_T_usage_list']

# Linear regression
model = LinearRegression()
model.fit(x, y)
y_pred = model.predict(x)

# Calculate R2
r2 = r2_score(y, y_pred)
print(f"R2 Score: {r2:.2f}")









'''import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

# Example DataFrame: A and T usage across 5000 bacteria


# Scatter plot to visualize the relationship
plt.scatter(file['P_A_usage_list'], file['C_T_usage_list'], alpha=0.5,color = 'green',s = 0.3)
plt.scatter(file['P_A_usage_list'], file['C_A_usage_list'], alpha=0.5, color=  "0.4",s = 0.3)
plt.xlabel('Ap Usage')
plt.ylabel('Ac and Tc Nucleotide Usage')
plt.title('Scatter Plot of comparition Ac and Tc with Ap Usage')
plt.legend(['Tc', 'Ac'],markerscale=15, handlelength=3, handleheight=2, fontsize=12)

plt.show()




# Pearson correlation coefficient
pearson_corr, _ = pearsonr(file['P_A_usage_list'], file['C_T_usage_list'])
print(f"Pearson Correlation Coefficient: {pearson_corr:.2f}")

# Spearman correlation coefficient
spearman_corr, _ = spearmanr(file['P_A_usage_list'], file['C_T_usage_list'])
print(f"Spearman Correlation Coefficient: {spearman_corr:.2f}")'''


