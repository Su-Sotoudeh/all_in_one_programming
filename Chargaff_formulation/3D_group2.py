import numpy as np
import pandas as pd
import plotly.graph_objects as go

# -------------------------
# 1) Load full dataset
# -------------------------

file_path = r'C:\Genomics\Chargaff\coding\all_in_one_programming\whole_datafram\Whole_datafram_2026_Eu_Pro_Virus_Plasmide.csv'
df_raw0 = pd.read_csv(file_path, low_memory=False)
print(df_raw0.head().to_string())

# فیلتر اولیه
df_raw = df_raw0[df_raw0['P2C_W'] < 8].copy()
groups = df_raw['Group'].unique()
print("Groups found:", groups)

# -------------------------
# 2) Precompute the green volume ONCE
# -------------------------

grid = 30  # smaller grid = faster WebGL
P2C = np.linspace(0, 10, grid)
A2T = np.linspace(0, 2, grid)
G2C = np.linspace(0, 2, grid)

X, Y, Z = np.meshgrid(P2C, A2T, G2C, indexing='ij')

M1 = (X * Y + 1) / (X + Y)
M2 = (X * Z + 1) / (X + Z)
Mw = np.abs(M1 - 1) + np.abs(M2 - 1) + 1

# -------------------------
# 3) Loop over groups
# -------------------------

for group in groups:

    print(f"Processing group: {group}")

    df = df_raw[df_raw['Group'] == group].copy()

    # تبدیل امن به عدد
    P2C_real = pd.to_numeric(df["P2C_W"], errors='coerce').values
    A2T_real = pd.to_numeric(df["A2T_annot"], errors='coerce').values
    G2C_real = pd.to_numeric(df["G2C_annot"], errors='coerce').values

    # حذف NaN ها
    mask_valid = ~(
            np.isnan(P2C_real) |
            np.isnan(A2T_real) |
            np.isnan(G2C_real)
    )
    P2C_real = P2C_real[mask_valid]
    A2T_real = A2T_real[mask_valid]
    G2C_real = G2C_real[mask_valid]

    # محدود کردن به محدوده‌ی شبکه
    mask_range = (
            (P2C_real >= 0) & (P2C_real <= 10) &
            (A2T_real >= 0) & (A2T_real <= 2) &
            (G2C_real >= 0) & (G2C_real <= 2)
    )
    P2C_real = P2C_real[mask_range]
    A2T_real = A2T_real[mask_range]
    G2C_real = G2C_real[mask_range]

    # محاسبه‌ی Mw واقعی
    M1_real = (P2C_real * A2T_real + 1) / (P2C_real + A2T_real)
    M2_real = (P2C_real * G2C_real + 1) / (P2C_real + G2C_real)
    Mw_real = np.abs(M1_real - 1) + np.abs(M2_real - 1) + 1

    # انتخاب نقاط داخل محدوده‌ی سبز
    mask_real = (Mw_real >= 0.5) & (Mw_real <= 2)
    P2C_red = P2C_real[mask_real]
    A2T_red = A2T_real[mask_real]
    G2C_red = G2C_real[mask_real]

    # اگر نقطه‌ای نبود، از این گروه رد شو
    if len(P2C_red) == 0:
        print(f"No valid points in green region for group: {group}")
        continue

    # -------------------------
    # 4) Plot for this group
    # -------------------------

    fig = go.Figure()

    # حجم سبز
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

    # نقاط این گروه
    fig.add_trace(go.Scatter3d(
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
    ))

    # ثابت نگه داشتن محور‌ها
    fig.update_layout(
        scene=dict(
            xaxis=dict(title="+/-", range=[0, 10]),
            yaxis=dict(title="A/T", range=[0, 2]),
            zaxis=dict(title="G/C", range=[0, 2])
        ),
        title=f"Group: {group} (Mw = 0.5–2)"
    )

    # اگر خواستی نقاط روی حجم باشند:
    fig.data = fig.data[::-1]

    fig.show()
