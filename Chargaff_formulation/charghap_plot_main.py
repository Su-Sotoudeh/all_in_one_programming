'''Ap = 10/100#0.1 /(0.1+0.4+0.8+0.2)
Tp = 20/100#0.4 /(0.1+0.4+0.8+0.2)
Ac = 46.66/100#0.8 /(0.1+0.4+0.8+0.2)
Tc = 23.33/100#0.2 /(0.1+0.4+0.8+0.2)

x = (Ap+Tp)/(Ac+Tc)
x1 = (Ap+Tp)/(Ac+Tc)
x2 = Tp/Ac
x3 = Ap/Tc
print('P2C',x1,x2,x3)
y = Ap/Tp
print('A2T',y)
M = (Ap+Ac)/(Tp+Tc)
#M2 =((1/(x+1)*(y+1))+ (y/(y+1))-(y/(x+1)*(y+1)))/((y/(x+1)*(y+1))+(1/(y+1))-(1/(x+1)*(y+1)))
M3 = ((y*x) +1) / (x+y)
print('M',M)
#print('M2',M2)
print('M3',M3)'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


'''file = pd.read_excel(r'D:\PhD\chargaf_genome_test_project\bacteria\CDS_2024-12-29_00-14-30\results_2024-12-29_00-14-30.xlsx')




# Define the range for x and y
x = np.linspace(0.0, 3, 1000)
y = np.linspace(0.0, 3, 1000)

# Create a meshgrid for x and y
X, Y = np.meshgrid(x, y)

# Define the new formula for M
M = ((Y * X) + 1) / (X + Y)

# User-provided lists for x and y values
list_x = file['P2C']
list_y = file['A2T_P']

# Calculate M values for the provided lists
list_M = [((y * x) + 1) / (x + y) for x, y in zip(list_x, list_y)]

# Plot the heatmap for M
plt.figure(figsize=(10, 8))  # Adjust the figure size
cmap = plt.cm.viridis  # Use a perceptually uniform colormap
contour = plt.contourf(X, Y, M, levels=np.linspace(0.7, 1.3, 50), cmap=cmap)

# Plot the user-provided (x, y) points on the heatmap as black dots
plt.scatter(list_x, list_y, color='black', label="Provided Data Points", edgecolor='white')

# Add colorbar and labels
cbar = plt.colorbar(contour, pad=0.02, aspect=30)  # Adjust colorbar placement and size
cbar.set_label("M Value")
plt.title("Heatmap of M with Provided Data Points")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.show()'''




import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap
import pandas as pd

#file = pd.read_excel(r'D:\PhD\chargaf_genome_test_project\bacteria\CDS_2024-12-29_00-14-30\results_2024-12-29_00-14-30.xlsx')
file = pd.read_excel(r'C:\Genomics\Chargaff\coding\all_in_one_programming\ML\ML_Bac\CDS_2024-12-29_00-14-30\results_2024-12-29_00-14-30.xlsx')


#chordata = pd.read_excel(r'D:\PhD\chargaf_genome_test_project\chordata\refrence_CDS_2025-01-13_06-29-06\results_2025-01-13_06-29-06.xlsx')


#plant = pd.read_excel(r"D:\Projects\all_in_one_programming\some_result\plant\results_2025-01-20_00-29-52.xlsx")

#insect = pd.read_excel(r"D:\Projects\all_in_one_programming\some_result\insect\refrence_CDS_2025-02-02_05-57-17\results_2025-02-02_05-57-17.xlsx")
#insect = insect[(insect["A2T_W"]<0.95) | (insect['A2T_W']>1.05)]
#print(insect.to_string())
# Define the range for x and y
x = np.linspace(0, 20, 1000)
y = np.linspace(0, 5, 1000)

# Create a meshgrid for x and y
X, Y = np.meshgrid(x, y)

# Define the new formula for M
M = ((Y * X) + 1) / (X + Y)

# User-provided lists for x and y values
list_x = file['P2C']
list_y = file['A2T_P']
'''chordatax = chordata['P2C']
chordatay = chordata['A2T_P']


plantx = plant["P2C"]
planty = plant["A2T_P"]
insectx = insect["P2C"]
insecty = insect["A2T_P"]'''
#insectx = ["1","1.5","2"]
#insecty = ["1","1.5","2"]



# Create a custom colormap with green near 1 and red farther away
colors = [(1, 0, 0), (0, 1, 0), (1, 0, 0)]  # Red, Green, Red
cmap = LinearSegmentedColormap.from_list("custom_red_green", colors)

# Set normalization with the desired range
norm = Normalize(vmin=0.95, vmax=1.05)        #green region

# Plot the heatmap for M
plt.figure(figsize=(6, 6))  # Adjust the figure size
contour = plt.contourf(X, Y, M, levels=np.linspace(0.90, 1.1, 1000), cmap=cmap, norm=norm)#red region

# Plot the user-provided (x, y) points on the heatmap as black dots
# Add contour lines to indicate M values
# Add contour lines to indicate specific M values
contour_levels = [0.90, 0.95, 1.0, 1.05, 1.1]  # Customize the levels as needed
contour_lines = plt.contour(X, Y, M, levels=contour_levels, colors='black', linewidths=0.7)

# Label the contour lines with their values
plt.clabel(contour_lines, inline=True, fontsize=8, fmt="M=%.2f")


#plt.scatter(chordatax, chordatay, color='w',s = 10,alpha=0.5, label="Provided Data Points", edgecolor='black')

# Add colorbar and labels
cbar = plt.colorbar(contour, pad=0.02, aspect=30)
cbar.set_label("M Value")
#cbar.set_clim(0.99, 1.01)  # Force colorbar limits

plt.title("Heatmap of M with Provided Data Points")
plt.xlabel("P2C")
plt.ylabel("A2T")
plt.legend()
plt.show()



