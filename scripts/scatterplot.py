#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, zscore
from sklearn.linear_model import LinearRegression
import pandas as pd

df = pd.read_csv("databases/phychem.tsv", sep = "\t")
df.sort_values(by='Molecular_Weight', ascending=False).head()

# CONFIG
X_COL = "Isoelectric_Point"
Y_COL = "Amphipathicity"

OUTLIER_Z_THRESHOLD = 1.5

# PREPARE DATA
plot_df = df[
    [X_COL, Y_COL, "Net_Charge", "Molecular_Weight", "ID", "Sequence"]
].copy()

plot_df = plot_df.dropna(subset=[X_COL, Y_COL])

# PEARSON CORRELATION
r, p_value = pearsonr(plot_df[X_COL], plot_df[Y_COL])

# LINEAR REGRESSION + RESIDUALS
X = plot_df[[X_COL]].values
y = plot_df[Y_COL].values

model = LinearRegression()
model.fit(X, y)

plot_df["Y_PRED"] = model.predict(X)
plot_df["RESIDUAL"] = plot_df[Y_COL] - plot_df["Y_PRED"]
plot_df["RESIDUAL_Z"] = zscore(plot_df["RESIDUAL"])

# OUTLIER CUTOFF
plot_df["OUTLIER"] = plot_df["RESIDUAL_Z"].abs() >= OUTLIER_Z_THRESHOLD

# PLOT
plt.rcParams.update({"font.size": 20})
sns.set_style("white")

plt.figure(figsize=(12, 8))

normal_df = plot_df[~plot_df["OUTLIER"]]
outlier_df = plot_df[plot_df["OUTLIER"]]

sns.scatterplot(
    data=normal_df,
    x=X_COL,
    y=Y_COL,
    hue="Net_Charge",
    size="Molecular_Weight",
    palette="coolwarm",
    sizes=(50, 300),
    alpha=0.65,
    legend=True
)

sns.scatterplot(
    data=outlier_df,
    x=X_COL,
    y=Y_COL,
    color="black",
    s=220,
    marker="X",
    label=f"Outlier |z| >= {OUTLIER_Z_THRESHOLD}",
    legend=True
)

sns.regplot(
    data=plot_df,
    x=X_COL,
    y=Y_COL,
    scatter=False,
    color="black",
    line_kws={"linewidth": 2, "linestyle": "--"}
)

plt.xlabel("Isoelectric Point", fontsize=20)
plt.ylabel("Amphipathicity", fontsize=20)

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

p_text = f"p = {p_value:.3f}" if p_value >= 0.001 else "p < 0.001"

plt.text(
    0.05,
    0.95,
    f"r = {r:.2f}\n{p_text}",
    transform=plt.gca().transAxes,
    fontsize=20,
    verticalalignment="top",
    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", lw=1)
)

plt.legend(
    bbox_to_anchor=(1.05, 1),
    loc=2,
    borderaxespad=0.,
    fontsize=16,
    title_fontsize=16
)

sns.despine()
plt.tight_layout()

plt.savefig(
    "results/scatterplot.svg",
    dpi=600,
    format="svg",
    bbox_inches="tight"
)

plt.close()
