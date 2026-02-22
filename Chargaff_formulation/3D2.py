import numpy as np
import pandas as pd
import plotly.graph_objects as go

# -------------------------
# 1) Load real dataset
# -------------------------

file_path = r"C:\Genomics\results\Genes_chardaff_parameters\Bacterial_2025-09-23_20-16-24.csv"

df = pd.read_csv(file_path)

# extract required columns
P2C_real = df["CDS_P2C"].values
A2T_real = df["CDS_P_A2T"].values
G2C_real = df["CDS_P_G2C"].values

# -------------------------
# 2) Compute M1, M2, Mw for real data
# -------------------------

M1_real = (P2C_real * A2T_real + 1) / (P2C_real + A2T_real)
M2_real = (P2C_real * G2C_real + 1) / (P2C_real + G2C_real)
Mw_real = np.abs(M1_real - 1) + np.abs(M2_real - 1) + 1

# select ONLY points inside lower green region (0.9–1.0)
mask_real_red = (Mw_real >= 0.95) & (Mw_real <= 1.05)

P2C_red = P2C_real[mask_real_red]
A2T_red = A2T_real[mask_real_red]
G2C_red = G2C_real[mask_real_red]

# -------------------------
# 3) Create the green volume (same as before)
# -------------------------

grid = 30
P2C = np.linspace(0, 3, grid)
A2T = np.linspace(0, 2, grid)
G2C = np.linspace(0, 2, grid)

X, Y, Z = np.meshgrid(P2C, A2T, G2C, indexing='ij')

M1 = (X * Y + 1) / (X + Y)
M2 = (X * Z + 1) / (X + Z)
Mw = np.abs(M1 - 1) + np.abs(M2 - 1) + 1

# -------------------------
# 4) Plot volume + real red points
# -------------------------

fig = go.Figure()

# green volume
fig.add_trace(go.Volume(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z.flatten(),
    value=Mw.flatten(),
    isomin=0.95,
    isomax=1.05,
    opacity=0.10,
    surface_count=40,
    colorscale=[[0, "green"], [1, "green"]]
))

# real test points (red)
fig.add_trace(go.Scatter3d(
    x=P2C_red,
    y=A2T_red,
    z=G2C_red,
    mode='markers',
    marker=dict(
        size=1,
        color='red',
        opacity=0.5
    ),
    name="Real Data Points"
))

fig.update_layout(
    scene=dict(
        xaxis_title="+/-",
        yaxis_title="A/T",
        zaxis_title="G/C"
    ),
    title="Real Genomic Data Inside Green Volume (Mw = 0.95–1.05)"
)

fig.show()
