import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Load your data
file_path = r'D:\PhD\chargaf_genome_test_project\bacteria\CDS_2024-12-28_23-01-19\results_2024-12-28_23-01-19.xlsx'
file_path2 = r'D:\PhD\chargaf_genome_test_project\chordata\refrence_CDS_2025-01-13_06-29-06\results_2025-01-13_06-29-06.xlsx'
file_path3 = r'D:\PhD\chargaf_genome_test_project\fungi2_git\fungi\refrence_CDS_2025-01-08_12-32-10\results_2025-01-08_12-32-10.xlsx'

path_list = [file_path,file_path2]
def A2T(file_path):
    file = pd.read_excel(file_path)

    # Correct way to select columns
    x = file.iloc[:, 1]
    y1 = file.iloc[:, 23]
    y2 = file.iloc[:, 29]
    y3 = file.iloc[:,14]
    y4 = file.iloc[:,5]

    # Calculate the range of x
    x_min = x.min()
    x_max = x.max()

    # Generate 10% of the total x values as evenly spaced points
    num_points = int(0.05 * len(x))  # 10% of total data
    x_points = np.linspace(x_min, x_max, num_points)


# Find the closest x values in the dataset to the evenly spaced points
    selected_indices = []
    for point in x_points:
        closest_idx = (np.abs(x - point)).argmin()
        selected_indices.append(closest_idx)

    # Create the scatter plot with smaller dots (using s parameter)
    plt.scatter(x, y1, color='blue', label='A2T_P', s=1)  # s=10 makes the dots smaller
    plt.scatter(x, y2, color='red', label='A2T_C', s=1)
    plt.scatter(x, y3, color='Green', label='A2T_PC', s=1)
    plt.scatter(x, y4, color='black', label='A2T_Wh', s=0.2)

    # Draw lines for the selected indices with evenly distributed points
    for i in selected_indices:
        plt.plot([x.iloc[i], x.iloc[i]], [y1.iloc[i], y2.iloc[i]], color='0.2', alpha=0.5, linewidth=1)

    # Add labels and legend
    plt.xlabel('A2T_W')
    plt.ylabel('A2T P,C,PC')
    plt.legend()
    plt.title('bacteria A2T')

    plt.show()




A2T(file_path2)