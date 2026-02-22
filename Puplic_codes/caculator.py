#x = P2C
#y = A2T_P
from fontTools.misc.arrayTools import pointInRect

'''x = 10
y = 1
A  =( (x/(x+1))  *  (y/y+1)        )   +       (    (1-(x/(x+1))) * (1-(y/(y+1)))    )
T  =( (x/(x+1))  *  (1-(y/(y+1)))  )   +       (    (1-(x/(x+1))) *  (y/y+1)         )

chargaf = A/T
print(chargaf)'''
#CDS_P2C 0.9107 CDS_P_A2T 1.1366 A2T_W 0.9984 Length(M) 3.623569
P2C = 0.9107
A2T = 1.1366


M = ((A2T * P2C) + 1) / (A2T + P2C)
#print(M)



'''A  =( (P2C/(P2C+1))  *  (A2T/(A2T+1))        )   +       (    (1-(P2C/(P2C+1))) * (1-(A2T/(A2T+1)))    )
print(A)
T  =( (P2C/(P2C+1))  *  (1-(A2T/(A2T+1)))  )   +       (    (1-(P2C/(P2C+1))) *  (A2T/(A2T+1))         )
print(T)
print(A/T)

X = 50000
Y = 1
A  =( (X/(X+1))  *  (Y/(Y+1))        )   +       (    (1-(X/(X+1))) * (1-(Y/(Y+1)))    )
T  =( (X/(X+1))  *  (1-(Y/(Y+1)))  )   +       (    (1-(X/(X+1))) *  (Y/(Y+1))         )
M = A/T
print(M)'''



import pandas as pd
import joblib

# -----------------------------
# LOAD MODEL + SCALER
# -----------------------------
model_path = r"C:\Genomics\Chargaff\ML_model_A2T.pkl"
scaler_path = r"C:\Genomics\Chargaff\y_scaler_A2T.pkl"

rf_model = joblib.load(model_path)
y_scaler = joblib.load(scaler_path)

# -----------------------------
# LOAD YOUR CSV DATA
# -----------------------------
df = pd.read_csv(r"C:\Genomics\results\Genes_chardaff_parameters\Whole_csv_file.csv")

A2T_arr = df["A2T_P"].values        # A2T
P2C_arr = df["P2C_W"].values        # P2C
Length_arr = df["Length(M)"].values # Length(M)
y_real = df["A2T_annot"].values     # real A2T annotation

# -----------------------------
# LOOP THROUGH ROWS
# -----------------------------
real = []
formula = []
ML = []
counter = 0
for P2C, A2T, L, real_y in zip(P2C_arr, A2T_arr, Length_arr, y_real):
    counter += 1
    if counter == 1000:
        break

    new_x = pd.DataFrame({
        "P2C_w": [P2C],
        "P_A2T_w": [A2T],
        "Length_M": [L]
    })

    # ML Prediction (scaled → original)
    pred_scaled = rf_model.predict(new_x.values)
    pred_original = y_scaler.inverse_transform(pred_scaled.reshape(-1, 1)).ravel()[0]

    # Mathematical M calculation
    M = ((A2T * P2C) + 1) / (A2T + P2C)

    # -----------------------------
    # PRINT RESULTS (ONLY FLOATS)
    # -----------------------------

    '''print(
        float(P2C),
        float(A2T),
        float(L),
        float(real_y),
        float(pred_original),
        float(M)
    )'''
    real.append(real_y)
    formula.append(pred_original)
    ML.append(M)
    print(round(real_y,3),round(pred_original,3),round(M,3))

import matplotlib.pyplot as plt
plt.plot(ML)
plt.plot(real)

plt.plot(real)
plt.show()










