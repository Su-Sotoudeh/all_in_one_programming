import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap
import pandas as pd

file = pd.read_excel(r'C:\Genomics\Chargaff\coding\all_in_one_programming\ML\ML_Bac\CDS_2024-12-29_00-14-30\results_2024-12-29_00-14-30.xlsx')


#file2 = pd.read_excel(r'D:\PhD\chargaf_genome_test_project\chordata\refrence_CDS_2025-01-13_06-29-06\results_2025-01-13_06-29-06.xlsx')

# Define the range for x and y
x = np.linspace(0, 2.5, 1000)
y = np.linspace(0, 2.5, 1000)

# Create a meshgrid for x and y
X, Y = np.meshgrid(x, y)

# Define the new formula for M
A  =( (X/(X+1))  *  (Y/(Y+1))        )   +       (    (1-(X/(X+1))) * (1-(Y/(Y+1)))    )
T  =( (X/(X+1))  *  (1-(Y/(Y+1)))  )   +       (    (1-(X/(X+1))) *  (Y/(Y+1))         )
M = A/T
print()
# User-provided lists for x and y values
list_x = file['P2C']
list_y = file['A2T_P']
'''list_x2 = file2['P2C']
list_y2 = file2['A2T_P']
list_x3 = file2['P2C']
list_y3 = file2['A2T_C']
list_x4 = file2['P2C']
list_y4 = file2['C2G_P']
list_x5 = file2['P2C']
list_y5 = file2['C2G_C']
list_x6 = file2['P2C']
list_y6 = file2['A2T_W']'''



# Calculate M values for the provided lists
#list_M = [((y * x) + 1) / (x + y) for x, y in zip(list_x, list_y)]

# Create a custom colormap with green near 1 and red farther away
colors = [(1, 0, 0), (0, 1, 0), (1, 0, 0)]  # Red, Green, Red
cmap = LinearSegmentedColormap.from_list("custom_red_green", colors)

# Set normalization with the desired range
norm = Normalize(vmin=0.90, vmax=1.1)

# Plot the heatmap for M
plt.figure(figsize=(6, 6))  # Adjust the figure size
contour = plt.contourf(X, Y, M, levels=np.linspace(0.80, 1.2, 1000), cmap=cmap, norm=norm)

# Plot the user-provided (x, y) points on the heatmap as black dots
# Add contour lines to indicate M values
# Add contour lines to indicate specific M values
contour_levels = [0.8, 0.95, 1.0, 1.05, 1.2]  # Customize the levels as needed
contour_lines = plt.contour(X, Y, M, levels=contour_levels, colors='black', linewidths=0.7)

# Label the contour lines with their values
plt.clabel(contour_lines, inline=True, fontsize=8, fmt="M=%.2f")

plt.scatter(list_x, list_y, color='red',s = 10,alpha=0.5, label="Provided Data Points", edgecolor='black')
#plt.scatter(list_x2, list_y2, color='blue',s = 10,alpha=0.5, label="Provided Data Points", edgecolor='blue')
#plt.scatter(list_x3, list_y3, color='yellow',s = 10,alpha=0.5, label="Provided Data Points", edgecolor='yellow')
#plt.scatter(list_x4, list_y4, color='red',s = 10,alpha=0.5, label="Provided Data Points", edgecolor='red')
#plt.scatter(list_x5, list_y5, color='pink',s = 10,alpha=0.5, label="Provided Data Points", edgecolor='pink')
#plt.scatter(list_x6, list_y2, color='black',s = 10,alpha=0.5, label="Provided Data Points", edgecolor='black')

# Add colorbar and labels
cbar = plt.colorbar(contour, pad=0.02, aspect=30)
cbar.set_label("M Value")
#cbar.set_clim(0.99, 1.01)  # Force colorbar limits

plt.title("Heatmap of M with Provided Data Points")
plt.xlabel("P2C")
plt.ylabel("A2T")
plt.legend()
plt.show()