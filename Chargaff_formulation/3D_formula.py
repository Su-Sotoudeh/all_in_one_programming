import numpy as np
import plotly.graph_objects as go

# 1) Define 2D axes
'''P2C = np.linspace(0, 2, 300)   # X axis
A2T = np.linspace(0, 2, 300)    # Y axis
X, Y = np.meshgrid(P2C, A2T)

# 2) Compute M1 = color map (2D)
M1 = (X * Y + 1) / (X + Y)

# 3) Define G2C axis (height levels)
G2C = np.linspace(0, 2, 5)   # 40 layers of the same 2D sheet

# 4) Colorscale: green visible, red transparent
transparent_red_green = [
    [0.0, "rgba(255,0,0,0)"],   # red → fully transparent
    [0.45, "rgba(255,0,0,0)"],  # red → fully transparent
    [0.49, "green"],            # green visible
    [0.51, "green"],            # green visible
    [0.55, "rgba(255,0,0,0)"],  # red → fully transparent
    [1.0, "rgba(255,0,0,0)"]    # red → fully transparent
]

# 5) Extrude the 2D sheet along G2C
surfaces = []
for g in G2C:
    surfaces.append(go.Surface(
        x=X,
        y=Y,
        z=np.ones_like(X) * g,     # same 2D sheet at height g
        surfacecolor=M1,           # SAME COLORS
        colorscale=transparent_red_green,
        cmin=0.90,
        cmax=1.10,
        showscale=False
    ))

fig = go.Figure(data=surfaces)

fig.update_layout(
    scene=dict(
        xaxis_title="P2C",
        yaxis_title="A2T",
        zaxis_title="G2C"
    ),
    title="Extruded 2D Heatmap Along G2C (Green Visible, Red Transparent)"
)

fig.show()'''





import numpy as np
import matplotlib.pyplot as plt

# 1) Generate datasets
N = 10000
A2T = np.random.uniform(0, 2, N)
G2C = np.random.uniform(0, 2, N)
P2C = np.random.uniform(0, 2, N)

# 2) Compute M1 = f(P2C, A2T)
M1 = (P2C * A2T + 1) / (P2C + A2T)

# 3) Define color rule
colors = np.where((M1 >= 0.95) & (M1 <= 1.05), "green", "red")

# 4) Plot 2D scatter (P2C vs A2T)
plt.figure(figsize=(7, 6))
plt.scatter(P2C, A2T, c=colors, s=5, alpha=0.7)

plt.xlabel("P2C")
plt.ylabel("A2T")
plt.title("2D M1 Map (Green = 0.95–1.05, Red = Outside Range)")

plt.show()



import numpy as np
import matplotlib.pyplot as plt

# 1) Generate datasets
N = 10000
A2T = np.random.uniform(0, 2, N)
P2C = np.random.uniform(0, 2, N)

# 2) Compute M1 = f(P2C, A2T)
M1 = (P2C * A2T + 1) / (P2C + A2T)

# 3) Find green points (M1 in range)
green_mask = (M1 >= 0.95) & (M1 <= 1.05)

P2C_green = P2C[green_mask]
A2T_green = A2T[green_mask]

# 4) Define possible G2C values
G2C_grid = np.linspace(0, 2, 200)

# 5) Build all combinations (vectorized)
P2C_rep = np.repeat(P2C_green, len(G2C_grid))
A2T_rep = np.repeat(A2T_green, len(G2C_grid))
G2C_rep = np.tile(G2C_grid, len(P2C_green))

# 6) Compute M2 for all combinations
M2 = (P2C_rep * G2C_rep + 1) / (P2C_rep + G2C_rep)

# 7) Keep only valid M2
valid_mask = (M2 >= 0.95) & (M2 <= 1.05)

P2C_3d = P2C_rep[valid_mask]
A2T_3d = A2T_rep[valid_mask]
G2C_3d = G2C_rep[valid_mask]

# 8) Plot 3D green points
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

ax.scatter(P2C_3d, A2T_3d, G2C_3d, c='green', s=5)

ax.set_xlabel("P2C")
ax.set_ylabel("A2T")
ax.set_zlabel("G2C")
ax.set_title("3D Green Points: M1 & M2 Both in Range")

plt.show()




