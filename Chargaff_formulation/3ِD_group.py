import numpy as np
import pandas as pd
import plotly.graph_objects as go

# -------------------------
# 1) Load full dataset
# -------------------------

file_path = r'C:\Genomics\Chargaff\coding\all_in_one_programming\whole_datafram\Whole_datafram_2026_Eu_Pro_Virus_Plasmide.csv'
df_raw0 = pd.read_csv(file_path, low_memory=False)
print(df_raw0.head().to_string())
df_raw = df_raw0[df_raw0['P2C_W']<8].copy()
groups = df_raw['Group'].unique()
print("Groups found:", groups)

# -------------------------
# 2) Precompute the green volume ONCE
# -------------------------

grid = 30  # smaller grid = faster WebGL
P2C = np.linspace(0, 10, grid)
A2T = np.linspace(0.8, 1.2, grid)
G2C = np.linspace(0.8, 1.2, grid)

X, Y, Z = np.meshgrid(P2C, A2T, G2C, indexing='ij')

M1 = (X * Y + 1) / (X + Y)
M2 = (X * Z + 1) / (X + Z)
Mw = np.abs(M1 - 1) + np.abs(M2 - 1) + 1

# -------------------------
# 3) Loop over groups
# -------------------------
counter = 0
for group in groups:
    if counter == 2:
        break

    print(f"Processing group: {group}")

    df =  df_raw[df_raw['Group'] == group]
    #df =  df_raw[(df_raw['Group'] != 'Plasmides')&(df_raw['Group'] != 'Viruses')&(df_raw['Group'] != 'Mito')&(df_raw['Group'] != 'Chloro_plast')]

    # extract columns
    P2C_real = df["P2C_W"].values
    A2T_real = df["A2T_annot"].values
    G2C_real = df["G2C_annot"].values


    # extract columns
    P2C_real = pd.to_numeric(df["P2C_W"], errors='coerce').values
    A2T_real = pd.to_numeric(df["A2T_annot"], errors='coerce').values
    G2C_real = pd.to_numeric(df["G2C_annot"], errors='coerce').values

    # Optional but very helpful: see where the problems are
    if df["P2C_W"].dtype == object or df["A2T_annot"].dtype == object:
        print(f"Problematic values in group {group}:")
        print("P2C_W non-numeric →", df["P2C_W"][~df["P2C_W"].astype(str).str.replace('.','',1).str.isnumeric()])
        print("A2T_annot non-numeric →", df["A2T_annot"][~df["A2T_annot"].astype(str).str.replace('.','',1).str.isnumeric()])




    # compute Mw
    M1_real = (P2C_real * A2T_real + 1) / (P2C_real + A2T_real)
    M2_real = (P2C_real * G2C_real + 1) / (P2C_real + G2C_real)
    Mw_real = np.abs(M1_real - 1) + np.abs(M2_real - 1) + 1

    # select points inside green region
    mask_real = (Mw_real >= 0.5) & (Mw_real <= 2)

    P2C_red = P2C_real[mask_real]
    A2T_red = A2T_real[mask_real]
    G2C_red = G2C_real[mask_real]






    # optional: limit number of points to avoid WebGL crash
    '''if len(P2C_red) > 3000:
        idx = np.random.choice(len(P2C_red), 3000, replace=False)
        P2C_red = P2C_red[idx]
        A2T_red = A2T_red[idx]
        G2C_red = G2C_red[idx]'''

    # -------------------------
    # 4) Plot for this group
    # -------------------------

    fig = go.Figure()

    # green volume
    fig.add_trace(go.Volume(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=Mw.flatten(),
        isomin=0.98,
        isomax=1.02,
        opacity=0.15,
        surface_count=20,
        colorscale=[[0, "green"], [1, "green"]],
        name="Green Volume"
    ))

    # red points for this group
    '''fig.add_trace(go.Scatter3d(
        x=P2C_red,
        y=A2T_red,
        z=G2C_red,
        mode='markers',
        marker=dict(
            size=1,
            color='black',
            opacity=0.7
        ),
        name=f"{group} points"
    ))'''


    fig.update_layout(
        scene=dict(
            xaxis_title="+/-",
            yaxis_title="A/T",
            zaxis_title="G/C"
        ),
        title=f"Group: {group} (Mw = 0.95–1.05)"
    )

    fig.show()
