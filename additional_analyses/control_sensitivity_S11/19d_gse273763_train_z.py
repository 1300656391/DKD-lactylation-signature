# 19d_gse273763_train_z.py ----
# Sensitivity companion for the GSE273763 FFPE ROC (Fig7A seventh curve):
# re-scores the submitter-normalized FFPE matrix with TRAINING-SET PDK2/PDK4
# means and SDs as the z anchor (instead of the within-cohort zrow used by
# 13_figures.R), keeping the frozen canonical coefficients. The within-z path
# must reproduce the pipeline's assert-backed AUC 0.805.
# Mirrors 13_figures.R L402-414 (matrix loading, labels, two-gene score) and
# 19e_export_train_mean_sd.R (training z anchor). Run from the project root.
# Output: results/revision/M5d_gse273763_train_z.csv
import csv
import gzip
import os

import numpy as np
from scipy.stats import mannwhitneyu

BETA = {"int": 1.85980920226601, "PDK2": -1.76986785255571, "PDK4": -2.43689946449653}
TRAIN_ANCHOR = {"PDK2": (5.952584, 0.437891), "PDK4": (9.863262, 1.838287)}  # from 19e (n=61)

src = r'data\GSE273763\GSE273763_DKD_GRCH38_forGEO_July2026_v2.txt.gz'
with gzip.open(src, "rt", encoding="utf-8", errors="replace") as f:
    import pandas as pd
    M = pd.read_csv(f, sep="\t", header=0, index_col=0)
M.columns = [str(c).replace(" ", "") for c in M.columns]
M = M.apply(lambda c: c.astype("float64") if c.dtype != "float64" else c)
assert M.shape[1] == 63, M.shape

y = np.array([1 if "DKD" in str(c) else 0 for c in M.columns], dtype=int)
assert y.sum() == 50 and (y == 0).sum() == 13


def auc(score):
    score = np.asarray(score, dtype=float)
    u = mannwhitneyu(score[y == 1], score[y == 0], alternative="two-sided")
    return u.statistic / (y.sum() * (y == 0).sum())


def zrow(v):  # R convention: sd with n-1 denominator
    v = np.asarray(v, dtype=float)
    return (v - v.mean()) / v.std(ddof=1)


def ztrain(v, gene):
    mu, sd = TRAIN_ANCHOR[gene]
    return (np.asarray(v, dtype=float) - mu) / sd


rows = []
for gene in ["PDK2", "PDK4"]:
    v = M.loc[gene].values.astype(float)
    a_w = auc(zrow(v))
    a_t = auc(ztrain(v, gene))
    rows.append(dict(feature=gene + " (single gene)", AUC_within_z=round(a_w, 3),
                     AUC_training_z=round(a_t, 3), delta=round(a_t - a_w, 3)))

s_w = (BETA["int"] + BETA["PDK2"] * zrow(M.loc["PDK2"].values.astype(float))
       + BETA["PDK4"] * zrow(M.loc["PDK4"].values.astype(float)))
s_t = (BETA["int"] + BETA["PDK2"] * ztrain(M.loc["PDK2"].values.astype(float), "PDK2")
       + BETA["PDK4"] * ztrain(M.loc["PDK4"].values.astype(float), "PDK4"))
auc_w = auc(s_w)
auc_t = auc(s_t)
rows.append(dict(feature="PDK2+PDK4 signature", AUC_within_z=round(auc_w, 3),
                 AUC_training_z=round(auc_t, 3), delta=round(auc_t - auc_w, 3)))

# pipeline-level assert: within-z must reproduce the shipped 0.805
assert abs(round(auc_w, 3) - 0.805) < 0.001, f"within-z AUC {auc_w:.4f} != 0.805"

out = r'results\revision\M5d_gse273763_train_z.csv'
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["feature", "AUC_within_z", "AUC_training_z", "delta"])
    w.writeheader()
    w.writerows(rows)
print(open(out, encoding="utf-8").read())
print("TRAINZ273_DONE")
