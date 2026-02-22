import pandas as pd
z = "/media/su/EC3A11743A113D50/genome_analysis/fungi/Bacterial_refrence_Int_UTR_CDS_2025-09-23_20-16-24/Bacterial_2025-09-23_20-16-24.xlsx"
c = pd.read_excel(z)
c.to_csv("/media/su/EC3A11743A113D50/genome_analysis/fungi/Bacterial_refrence_Int_UTR_CDS_2025-09-23_20-16-24/Bacterial_2025-09-23_20-16-24.csv",index=False)