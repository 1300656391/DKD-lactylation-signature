# evidence_gain/5_head_to_head/head_to_head.py ----
# 头对头比较：ref52（PTGS2/NFE2L2）vs 本文签名，同一批队列同一预处理。
# 口径说明：
#   - 子组分析（early/advanced）用子集内 z（与原管线一致，复得 0.889/1.000）
#   - 同时报告 as-trained 方向的原始 AUC 与 auto-oriented AUC（=max(A,1-A)，与 pROC direction='auto' 对齐，
#     即 Table 2 口径）；标记签名方向反转的语境
#   - 单基因按训练队列 rank-AUC 定向；mean-z 对按同定向合成；PTGS2/NFE2L2 冻结 logistic 镜像本文管线
import numpy as np, pandas as pd, os
from scipy.stats import mannwhitneyu

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "2_cross_deconv")
OUTD = os.path.dirname(os.path.abspath(__file__))
BETA = {"intercept": 1.85980920226601, "PDK2": -1.76986785255571, "PDK4": -2.43689946449653}
GENES = ["PDK2", "PDK4", "PTGS2", "NFE2L2"]

def load(coh):
    sub = pd.read_csv(os.path.join(BASE, f"bulk_{coh}_subset.csv"))
    ph = pd.read_csv(os.path.join(BASE, f"bulk_{coh}_pheno.csv"))
    expr = sub.set_index("gene"); expr = expr[~expr.index.duplicated(keep="first")]
    ph = ph.set_index("colname").loc[:, ["sample", "group"]]
    return expr[ph.index], ph

def zmat(e):
    return e.sub(e.mean(axis=1), axis=0).div(e.std(axis=1, ddof=1).replace(0, np.nan), axis=0)

def auc_raw(score, is_dn):
    u = mannwhitneyu(score[is_dn], score[~is_dn], alternative="two-sided")
    return u.statistic / (is_dn.sum() * (~is_dn).sum()), u.pvalue

def logistic_fit(X, y, iters=100):
    Xd = np.column_stack([np.ones(len(X)), X]); b = np.zeros(Xd.shape[1])
    for _ in range(iters):
        eta = np.clip(Xd @ b, -30, 30); p = 1 / (1 + np.exp(-eta))
        W = p * (1 - p) + 1e-9
        g = Xd.T @ (y - p); H = (Xd * W[:, None]).T @ Xd + 1e-8 * np.eye(Xd.shape[1])
        step = np.linalg.solve(H, g); b += step
        if np.max(np.abs(step)) < 1e-10: break
    return b

# ---- 训练队列：定向 + PTGS2/NFE2L2 冻结 logistic ----
e_tr, ph_tr = load("GSE96804")
Z_tr = zmat(e_tr.loc[GENES]); y_tr = (ph_tr["group"] != "Control").values.astype(int)
orient = {}
for g in GENES:
    a0, _ = auc_raw(Z_tr.loc[g].values, y_tr.astype(bool))
    orient[g] = 1 if a0 > 0.5 else -1
b_pair = logistic_fit(Z_tr.loc[["PTGS2", "NFE2L2"]].T.values, y_tr)
print("训练定向 orient:", orient)
print("PTGS2/NFE2L2 冻结 logistic: intercept=%.3f PTGS2=%.3f NFE2L2=%.3f" % tuple(b_pair))

CONTEXTS = [
    ("GSE96804_train",        "GSE96804",         None),
    ("GSE30528_glom",         "GSE30528",         None),
    ("GSE30529_tub",          "GSE30529",         None),
    ("GSE142025_allDN",       "GSE142025",        None),
    ("GSE142025_early",       "GSE142025",        "Early_DN"),
    ("GSE142025_advanced",    "GSE142025",        "Advanced_DN"),
    ("GSE104948_pooled",      "GSE104948_pooled", None),
    ("GSE104948_GPL22945",    "GSE104948_A",      None),
    ("GSE104948_GPL24120",    "GSE104948_B",      None),
]
TABLE2 = {"GSE30528_glom": 0.906, "GSE30529_tub": 0.650, "GSE142025_allDN": 0.992,
          "GSE142025_early": 0.889, "GSE142025_advanced": 1.000, "GSE104948_pooled": 1.000,
          "GSE104948_GPL22945": 1.000, "GSE104948_GPL24120": 0.600}

cache = {}; rows = []
for ctx, coh, subgrp in CONTEXTS:
    if coh not in cache:
        e, ph = load(coh); cache[coh] = (zmat(e.loc[GENES]), ph)
    Z, ph = cache[coh]
    if subgrp:  # 子集内 z（原管线口径）：先截取样本再 z
        keep = ((ph["group"] == subgrp) | (ph["group"] == "Control")).values
        sub = pd.read_csv(os.path.join(BASE, f"bulk_{coh}_subset.csv")).set_index("gene")
        sub = sub[~sub.index.duplicated(keep="first")]
        php = pd.read_csv(os.path.join(BASE, f"bulk_{coh}_pheno.csv")).set_index("colname").loc[:, ["sample", "group"]]
        sub = sub[php.index]
        Z = zmat(sub.loc[GENES].loc[:, keep]); ph = php[keep]
    is_dn = (ph["group"] != "Control").values.astype(bool)
    feats = {
        "PDK2":              orient["PDK2"] * Z.loc["PDK2"],
        "PDK4":              orient["PDK4"] * Z.loc["PDK4"],
        "PTGS2":             orient["PTGS2"] * Z.loc["PTGS2"],
        "NFE2L2":            orient["NFE2L2"] * Z.loc["NFE2L2"],
        "PTGS2+NFE2L2_meanz": (orient["PTGS2"] * Z.loc["PTGS2"] + orient["NFE2L2"] * Z.loc["NFE2L2"]) / 2,
        "PTGS2+NFE2L2_frozen": b_pair[0] + b_pair[1] * Z.loc["PTGS2"] + b_pair[2] * Z.loc["NFE2L2"],
        "PDK2/PDK4_signature": BETA["intercept"] + BETA["PDK2"] * Z.loc["PDK2"] + BETA["PDK4"] * Z.loc["PDK4"],
    }
    for name, score in feats.items():
        a, p = auc_raw(score.values.astype(float), is_dn)
        rows.append(dict(context=ctx, n_dn=int(is_dn.sum()), n_ctrl=int((~is_dn).sum()), feature=name,
                         AUC_as_trained=round(a, 3), AUC_oriented=round(max(a, 1 - a), 3),
                         direction_reversed=(a < 0.5), p=round(p, 6),
                         Table2=TABLE2.get(ctx) if name == "PDK2/PDK4_signature" else np.nan))
res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUTD, "head_to_head_auc.csv"), index=False)

# ---- 配对 bootstrap：签名 vs PTGS2（ref52 最强特征），B=2000，方向用 oriented 口径 ----
rng = np.random.default_rng(42)
bt_rows = []
for ctx, coh, subgrp in CONTEXTS:
    if ctx == "GSE96804_train":
        continue
    if coh not in cache:
        continue
    Z, ph = cache[coh]
    if subgrp:
        keep = ((ph["group"] == subgrp) | (ph["group"] == "Control")).values
        sub = pd.read_csv(os.path.join(BASE, f"bulk_{coh}_subset.csv")).set_index("gene")
        sub = sub[~sub.index.duplicated(keep="first")]
        php = pd.read_csv(os.path.join(BASE, f"bulk_{coh}_pheno.csv")).set_index("colname").loc[:, ["sample", "group"]]
        sub = sub[php.index]
        Z = zmat(sub.loc[GENES].loc[:, keep]); ph = php[keep]
    is_dn = (ph["group"] != "Control").values.astype(bool)
    s1 = (BETA["intercept"] + BETA["PDK2"] * Z.loc["PDK2"] + BETA["PDK4"] * Z.loc["PDK4"]).values
    s2 = (orient["PTGS2"] * Z.loc["PTGS2"]).values
    def or_auc(x):
        a, _ = auc_raw(x, is_dn); return max(a, 1 - a)
    a1, a2 = or_auc(s1), or_auc(s2)
    idx = np.arange(len(s1)); diffs = []
    for _ in range(2000):
        b = rng.choice(idx, size=len(idx), replace=True)
        bd = is_dn[b]
        if bd.sum() < 2 or (~bd).sum() < 2:
            continue
        ab = mannwhitneyu(s1[b][bd], s1[b][~bd]).statistic / (bd.sum() * (~bd).sum())
        cd = mannwhitneyu(s2[b][bd], s2[b][~bd]).statistic / (bd.sum() * (~bd).sum())
        diffs.append(max(ab, 1 - ab) - max(cd, 1 - cd))
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    p_two = 2 * min(np.mean(np.array(diffs) <= 0), np.mean(np.array(diffs) >= 0))
    bt_rows.append(dict(context=ctx, AUC_signature=round(a1, 3), AUC_PTGS2=round(a2, 3),
                        diff=round(a1 - a2, 3), boot_lo=round(lo, 3), boot_hi=round(hi, 3),
                        p_boot=round(min(p_two, 1.0), 4)))
bt = pd.DataFrame(bt_rows)
bt.to_csv(os.path.join(OUTD, "signature_vs_PTGS2_bootstrap.csv"), index=False)

print("\n== 本文签名复现检查（AUC_oriented 对照 Table 2） ==")
chk = res[res.feature == "PDK2/PDK4_signature"][["context", "AUC_as_trained", "AUC_oriented", "Table2", "direction_reversed"]]
print(chk.to_string(index=False))
print("\n== AUC 矩阵（AUC_oriented，Table 2 口径） ==")
print(res.pivot_table(index="feature", columns="context", values="AUC_oriented").round(3).to_string())
print("\n== 方向反转标记（as-trained AUC<0.5 的语境 x 特征） ==")
rev = res[res.direction_reversed][["context", "feature", "AUC_as_trained"]]
print(rev.to_string(index=False) if len(rev) else "（无）")
print("\n== 签名 vs PTGS2 配对 bootstrap（B=2000） ==")
print(bt.to_string(index=False))
print("\nHEADTOHEAD_DONE")
