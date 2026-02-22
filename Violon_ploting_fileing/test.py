import matplotlib
matplotlib.rcParams.update(matplotlib.rcParamsDefault)  # reset everything
matplotlib.rcParams['figure.dpi'] = 300                 # fixed DPI
matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['figure.figsize'] = (2,1)          # fixed figure size
matplotlib.rcParams['font.size'] = 10                  # base font
matplotlib.rcParams['axes.titlesize'] = 12
matplotlib.rcParams['axes.labelsize'] = 10
matplotlib.rcParams['xtick.labelsize'] = 8
matplotlib.rcParams['ytick.labelsize'] = 8
matplotlib.rcParams['lines.linewidth'] = 0.8
matplotlib.rcParams['axes.linewidth'] = 0.8

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib

# ---- FIX SCALING ----
matplotlib.rcParams.update(matplotlib.rcParamsDefault)
matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['figure.figsize'] = (14,2)
matplotlib.rcParams['font.size'] = 10
matplotlib.rcParams['axes.titlesize'] = 12
matplotlib.rcParams['axes.labelsize'] = 10
matplotlib.rcParams['xtick.labelsize'] = 8
matplotlib.rcParams['ytick.labelsize'] = 8
matplotlib.rcParams['lines.linewidth'] = 0.8
matplotlib.rcParams['axes.linewidth'] = 0.8

sns.reset_defaults()
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.0)

# ---- TEST DATA ----
df = pd.DataFrame(np.random.rand(10,4)*100, columns=["A","T","C","G"])

# ---- PLOT ----
fig, ax = plt.subplots()
sns.violinplot(data=df, ax=ax, inner="box", cut=0, linewidth=0.8)
ax.set_ylim(0,100)
plt.show()

