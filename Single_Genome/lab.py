import pandas as pd

file = r'D:\khalili.tsv'

df = pd.read_csv(file, sep='\t')
print(df.head())
