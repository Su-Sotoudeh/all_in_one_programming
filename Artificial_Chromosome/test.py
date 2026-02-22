import random
import numpy as np
import pandas as pd

print(random.randint(50,180))

a = pd.read_csv(r"C:\Genomics\results\Genes_chardaff_parameters\Bacteria_2025-11-11_17-47-44.csv")
A2T = a['CDS_P2C']