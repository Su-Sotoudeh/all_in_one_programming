import pandas as pd

file = pd.read_excel('RNA_plus.xlsx',sheet_name='Sheet1',usecols = "A")
datafram = pd.DataFrame(file)
print(datafram[1])
file = str(file)
file = file.split()
print(file)