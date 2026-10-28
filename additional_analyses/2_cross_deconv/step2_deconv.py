# evidence_gain/2_cross_deconv/step2_deconv.py ----
# 任务#2 交叉去卷积（与主稿 ssGSEA+Bindea 方法学独立）：
#   参照 = GSE131882 肾 snRNA 10 型线性 CP10K 均值（marker 矩阵，本步输入）
#   混合 = 各 bulk 队列表达（队列内 per-gene z 标准化，尺度无关，规避 log/线性与平台差异）
#   求解 = scipy NNLS（CIBERSORT 同类的参照约束去卷积；系数归一化为分数）
#   检验 = 各型分数 DN vs Ctrl（Wilcoxon, BH10）；PDK2/PDK4/签名分 vs 各型分数 Spearman（全样本 + DN 内）
# 注：GSE131882 参照不含 Macro/T-NK（0 核），B 细胞为唯一免疫型 → 免疫结论以 B 细胞分数为代理，
#     宏噬细胞特异复现不在此参照能力范围内（如实报告）。
import numpy as np, pandas as pd
from scipy.optimize import nnls
from scipy.stats import spearmanr, mannwhitneyu
from itertools import product
import os

BASE = os.path.dirname(os.path.abspath(__file__))
COHORTS = ["GSE96804", "GSE30528", "GSE30529", "GSE142025", "GSE104948_A", "GSE104948_B", "GSE104948_pooled"]
COHORT_LABEL = {"GSE104948_A": "GSE104948_GPL22945", "GSE104948_B": "GSE104948_GPL24120",
                "GSE104948_pooled": "GSE104948_pooled"}
BETA = {"intercept": 1.85980920226601, "PDK2": -1.76986785255571, "PDK4": -2.43689946449653}

ref = pd.read_csv(os.path.join(BASE, "ref_signature_marker_linear.csv"), index_col=0)
types = list(ref.columns)

def zscore(df):
    mu = df.mean(axis=1); sd = df.std(axis=1, ddof=0).replace(0, np.nan)
    return df.sub(mu, axis=0).div(sd, axis=0).dropna()

ref_z = zscore(ref)  # per-gene across types

all_frac, test_rows, corr_rows = [], [], []
for coh in COHORTS:
    sub = pd.read_csv(os.path.join(BASE, f"bulk_{coh}_subset.csv"))
    ph  = pd.read_csv(os.path.join(BASE, f"bulk_{coh}_pheno.csv"))
    expr = sub.set_index("gene")
    expr = expr[~expr.index.duplicated(keep="first")]
    ph = ph.set_index("colname").loc[:, ["sample", "group"]]
    expr = expr[ph.index]
    ez = zscore(expr)
    genes = ref_z.index.intersection(ez.index)
    R = ref_z.loc[genes, :].values
    fracs = {}
    for s in ez.columns:
        coef, _ = nnls(R, ez.loc[genes, s].values)
        fracs[s] = coef / coef.sum() if coef.sum() > 0 else np.full(len(types), np.nan)
    F = pd.DataFrame(fracs, index=types).T
    F["sample"], F["group"], F["cohort"] = ph["sample"].values, ph["group"].values, COHORT_LABEL.get(coh, coh)
    all_frac.append(F.reset_index(drop=True))

    is_dn = F["group"] != "Control"
    for t in types:
        x = F[t]; u = mannwhitneyu(x[is_dn], x[~is_dn], alternative="two-sided")
        test_rows.append(dict(cohort=COHORT_LABEL.get(coh, coh), type=t,
                               median_DN=float(x[is_dn].median()), median_Ctrl=float(x[~is_dn].median()),
                               direction="up" if x[is_dn].median() > x[~is_dn].median() else "down",
                               p=float(u.pvalue)))
    # PDK 表达与签名分（队列内 z）
    def gene_row(g): return expr.loc[g] if g in expr.index else None
    zp = {g: (gene_row(g) - gene_row(g).mean()) / gene_row(g).std(ddof=0) for g in ("PDK2", "PDK4") if gene_row(g) is not None}
    if len(zp) == 2:
        score = BETA["intercept"] + BETA["PDK2"] * zp["PDK2"] + BETA["PDK4"] * zp["PDK4"]
    else:
        score = None
    feats = {"PDK2": zp.get("PDK2"), "PDK4": zp.get("PDK4"), "signature_score": score}
    for fname, fv in feats.items():
        if fv is None: continue
        for t in types:
            r_all = spearmanr(fv, F[t])
            r_dn  = spearmanr(fv[is_dn.values], F[t][is_dn]) if is_dn.sum() >= 5 else (np.nan, np.nan)
            corr_rows.append(dict(cohort=COHORT_LABEL.get(coh, coh), feature=fname, type=t,
                                  rho_all=round(float(r_all.statistic), 3),
                                  p_all=float(r_all.pvalue),
                                  rho_within_DN=round(float(r_dn[0]), 3) if not np.isnan(r_dn[0]) else np.nan))

F = pd.concat(all_frac, ignore_index=True)
F.to_csv(os.path.join(BASE, "deconv_fractions_all.csv"), index=False)
tests = pd.DataFrame(test_rows)
for coh in tests["cohort"].unique():
    m = tests["cohort"] == coh
    tests.loc[m, "p_BH10"] = tests.loc[m, "p"].rank(method="average").mul(0)  # placeholder keep order
    n = int(m.sum())
    tests.loc[m, "p_BH10"] = tests.loc[m, "p"] * n / tests.loc[m, "p"].rank(method="average")
tests["p_BH10"] = tests["p_BH10"].clip(upper=1.0).round(6)
tests.to_csv(os.path.join(BASE, "deconv_group_tests.csv"), index=False)
corr = pd.DataFrame(corr_rows)
corr["p_all"] = corr["p_all"].round(6)
corr.to_csv(os.path.join(BASE, "deconv_pdk_correlations.csv"), index=False)

print("== fractions: median by cohort x group (B / PT / Podo) ==")
print(F.groupby(["cohort", "group"])[["B", "PT", "Podo"]].median().round(4).to_string())
print("\n== group tests (BH10 per cohort), p_BH10 < 0.05 ==")
print(tests[tests["p_BH10"] < 0.05][["cohort", "type", "median_DN", "median_Ctrl", "direction", "p", "p_BH10"]].round(4).to_string(index=False))
print("\n== PDK2 vs fractions (rho_all), primary cohorts ==")
prim = corr[(corr["feature"] == "PDK2") & (corr["cohort"].isin(["GSE96804", "GSE142025", "GSE30528", "GSE30529"]))]
print(prim[["cohort", "type", "rho_all", "p_all", "rho_within_DN"]].round(4).to_string(index=False))
print("\nDECONV_DONE")
