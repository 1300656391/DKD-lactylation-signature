# evidence_gain/4_new_cohort/score_gse273763_counts.py ----
# GSE273763 正式复评（审计第 10 条：归一化矩阵疑有 4095 封顶 -> 用 counts 矩阵过标准归一化）
# 主口径 log2(CPM+1)；敏感性 DESeq2 中位数比值法（median-of-ratios）+ log2
import gzip
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

BASE = r"C:\Users\cangf\ZCodeProject\糖肾\evidence_gain\4_new_cohort"
BETA = {"intercept": 1.85980920226601, "PDK2": -1.76986785255571, "PDK4": -2.43689946449653}

with gzip.open(BASE + r"\GSE273763_counts_symbol.csv.gz", "rt", encoding="utf-8", errors="replace") as f:
    df = pd.read_csv(f, sep=",", index_col=0)
df.index = df.index.astype(str).str.upper().str.strip()
cols = [str(c) for c in df.columns]
ph = pd.Series(["DKD" if "DKD" in c.upper() else "Donor" for c in cols], index=cols)
print("dim:", df.shape, "| groups:", dict(ph.value_counts()), "| 行名示例:", list(df.index[:3]))
if df.index.duplicated().any():
    mu = df.sum(axis=1)
    keep = df.index.to_series().groupby(df.index).apply(lambda i: i[np.argmax(mu.loc[i])])
    df = df.loc[keep.values]
    print("重复基因折叠后:", df.shape)
num = df.select_dtypes(include=[np.number])
print("整数计数比例: %.3f | max: %d | 4095 值计数: %d" % (
    (num.values == np.round(num.values)).mean(), num.values.max(), (num.values == 4095).sum()))

def norm_cpm(m):
    lib = m.sum(axis=0); lib[lib <= 0] = np.nan
    return np.log2(m.div(lib, axis=1) * 1e6 + 1)

def norm_mor(m):  # DESeq2 中位数比值
    lib = m.sum(axis=0); lib[lib <= 0] = np.nan
    g = m[m.sum(axis=1) > 0]
    ref = np.exp(np.log(g.div(lib, axis=1)).replace([np.inf, -np.inf], np.nan).median(axis=1))
    size = np.nanmedian(m.div(ref, axis=0), axis=0); size[size <= 0] = 1
    return np.log2(m.div(size, axis=1) + 1)

GENES = ["PDK2", "PDK4", "PTGS2", "NFE2L2"]
orient = {g: -1 for g in GENES}  # 训练队列定向（四基因在训练中均 DN 下调）
is_dn = (ph == "DKD").values

def auc_raw(score):
    u = mannwhitneyu(np.asarray(score)[is_dn], np.asarray(score)[~is_dn], alternative="two-sided")
    return u.statistic / (is_dn.sum() * (~is_dn).sum()), u.pvalue

all_rows = []
for tag, M in [("log2CPM", norm_cpm(num)), ("medianRatio", norm_mor(num))]:
    E = M.loc[[g for g in GENES if g in M.index]]
    Z = (E - E.mean(axis=1).values[:, None]) / E.std(axis=1, ddof=1).values[:, None]
    feats = {g: orient[g] * Z.loc[g] for g in Z.index}
    if {"PDK2", "PDK4"} <= set(feats):
        feats["PDK2/PDK4_signature"] = BETA["intercept"] + BETA["PDK2"] * Z.loc["PDK2"] + BETA["PDK4"] * Z.loc["PDK4"]
    if {"PTGS2", "NFE2L2"} <= set(feats):
        feats["PTGS2+NFE2L2_meanz"] = (feats["PTGS2"] + feats["NFE2L2"]) / 2
    for name, s in feats.items():
        a, p = auc_raw(s.values)
        all_rows.append(dict(norm=tag, feature=name, AUC_as_trained=round(a, 3),
                             AUC_oriented=round(max(a, 1 - a), 3), p=("%.3g" % p)))
    # 代理相关
    sig = feats.get("PDK2/PDK4_signature")
    for label, genes in {"tubular": ["LRP2", "CUBN"], "fibrosis": ["COL1A1", "COL3A1"]}.items():
        hit = [g for g in genes if g in M.index]
        if len(hit) == 2 and sig is not None:
            proxy = M.loc[hit].mean(axis=0)
            r, pv = spearmanr(sig.values, proxy.values)
            all_rows.append(dict(norm=tag, feature=f"signature x {label}", AUC_as_trained=round(r, 3),
                                 AUC_oriented=np.nan, p=("%.3g" % pv)))
out = pd.DataFrame(all_rows)
out.to_csv(BASE + r"\GSE273763_counts_scoring.csv", index=False)
print(out.to_string(index=False))
print("\nCOUNTS_REEVAL_DONE")
