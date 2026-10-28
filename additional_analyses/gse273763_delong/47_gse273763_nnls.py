# evidence_gain/4_new_cohort/step1b_immune_nnls_gate.py ----
# 执行①补全：相关 p 值 + NNLS B 代理（与 #2 同法：GSE131882 十型参照 + 422 marker）+ 预注册判定
import numpy as np, pandas as pd, os
from scipy.optimize import nnls
from scipy.stats import spearmanr, mannwhitneyu

BASE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(BASE, "..", "2_cross_deconv")

# ---- ssGSEA 相关 p 值（读 R 已导出的逐样本分数）----
rows = []
for tag in ("submitter", "log2cpm"):
    ss = pd.read_csv(os.path.join(BASE, f"gse273763_ssgsea_{tag}.csv"))
    bulk = pd.read_csv(os.path.join(BASE, f"gse273763_bulk_subset_{tag}.csv")).set_index("gene") \
        if os.path.exists(os.path.join(BASE, f"gse273763_bulk_subset_{tag}.csv")) else None
    # PDK 表达直接从矩阵取
    import gzip
    if tag == "submitter":
        with gzip.open(os.path.join(BASE, "GSE273763_DKD_GRCH38_forGEO_July2026_v2.txt.gz"), "rt", encoding="utf-8", errors="replace") as f:
            M = pd.read_csv(f, sep="\t", index_col=0)
    else:
        with gzip.open(os.path.join(BASE, "GSE273763_counts_symbol.csv.gz"), "rt", encoding="utf-8", errors="replace") as f:
            C = pd.read_csv(f, index_col=0)
        C.index = C.index.astype(str).str.upper()
        M = np.log2(C.div(C.sum(axis=0), axis=1) * 1e6 + 1)
    M.columns = M.columns.str.replace(r"\s+", "", regex=True)
    pdk2 = M.loc["PDK2"]; pdk2.index = M.columns
    ss = ss.set_index("sample").loc[M.columns]
    is_dn = (ss["group"] == "DKD").values
    for cell in ("Macrophages", "B_cells", "T_cells", "Th2_cells"):
        for label, mask in (("all", np.ones(len(ss), bool)), ("DN_only", is_dn)):
            r, p = spearmanr(pdk2.values[mask], ss[cell].values[mask])
            rows.append(dict(caliber=tag, cell=cell, scope=label, rho=round(r, 3), p=round(p, 6)))
cor = pd.DataFrame(rows)
cor.to_csv(os.path.join(BASE, "gse273763_pdk_corr_pvalues.csv"), index=False)
print("== PDK2 x ssGSEA cells（含 p）=="); print(cor.to_string(index=False))

# ---- NNLS B 代理（#2 同法）----
ref = pd.read_csv(os.path.join(REF, "ref_signature_marker_linear.csv"), index_col=0)
mk = list(ref.index)
def zscore(df):
    mu = df.mean(axis=1); sd = df.std(axis=1, ddof=0).replace(0, np.nan)
    return df.sub(mu, axis=0).div(sd, axis=0).dropna()
ref_z = zscore(ref)
nrows = []
for tag, Mm in (("submitter", M),):
    pass  # M 已在上方构造（log2cpm）；submitter 版本单独构造
with gzip.open(os.path.join(BASE, "GSE273763_DKD_GRCH38_forGEO_July2026_v2.txt.gz"), "rt", encoding="utf-8", errors="replace") as f:
    S = pd.read_csv(f, sep="\t", index_col=0)
S.columns = S.columns.str.replace(r"\s+", "", regex=True)
S.index = S.index.astype(str).str.upper()
S = S.loc[~S.index.duplicated(keep="first")]
grp = pd.Series(["DKD" if ("DKD" in c) else "Donor" for c in S.columns], index=S.columns)
for tag, E in (("submitter", S), ("log2cpm", M)):
    hit = [g for g in set(mk + ["PDK2", "PDK4"]) if g in E.index]
    sub = E.loc[hit]
    ez = zscore(sub)
    genes = ref_z.index.intersection(ez.index)
    R = ref_z.loc[genes, :].values
    fr = {}
    for s in ez.columns:
        coef, _ = nnls(R, ez.loc[genes, s].values)
        fr[s] = coef / coef.sum() if coef.sum() > 0 else np.full(ref.shape[1], np.nan)
    F = pd.DataFrame(fr, index=ref.columns).T
    F["group"] = grp.values
    F.to_csv(os.path.join(BASE, f"gse273763_nnls_fractions_{tag}.csv"))
    isd = (F["group"] == "DKD").values
    u = mannwhitneyu(F["B"][isd], F["B"][~isd], alternative="two-sided")
    r_all, p_all = spearmanr(sub.loc["PDK2"].values, F["B"].values)
    r_dn, p_dn = spearmanr(sub.loc["PDK2"].values[isd], F["B"].values[isd])
    nrows.append(dict(caliber=tag, B_median_DN=round(F["B"][isd].median(), 4), B_median_Ctrl=round(F["B"][~isd].median(), 4),
                      B_p=round(u.pvalue, 6), PDK2xB_rho=round(r_all, 3), PDK2xB_p=round(p_all, 6),
                      PDK2xB_DN_rho=round(r_dn, 3), PDK2xB_DN_p=round(p_dn, 6)))
nn = pd.DataFrame(nrows)
nn.to_csv(os.path.join(BASE, "gse273763_nnls_B_proxy.csv"), index=False)
print("\n== NNLS B 代理（#2 同法，十型参照）=="); print(nn.to_string(index=False))

# ---- 预注册判定（锚点 B 节；主口径=submitter）----
pri = cor[(cor.caliber == "submitter")]
mac_pri = pri[(pri.cell == "Macrophages") & (pri.scope == "all")].iloc[0]
mac_up = (mac_pri.rho < 0)
tests = pd.read_csv(os.path.join(BASE, "gse273763_immune_tests_submitter.csv"))
mac_bh = tests[tests.cell == "Macrophages"].iloc[0]
gate_primary = (mac_bh["direction"] == "up") and (mac_bh["p.adj_BH18"] < 0.05) and mac_up and (mac_pri.p < 0.05)
sens = cor[(cor.caliber == "log2cpm")]
mac_sens = sens[(sens.cell == "Macrophages") & (sens.scope == "all")].iloc[0]
tests_s = pd.read_csv(os.path.join(BASE, "gse273763_immune_tests_log2cpm.csv"))
mac_bh_s = tests_s[tests_s.cell == "Macrophages"].iloc[0]
gate_sens = (mac_bh_s["direction"] == "up") and (mac_bh_s["p.adj_BH18"] < 0.05) and (mac_sens.rho < 0) and (mac_sens.p < 0.05)
print("\n== 预注册判定 ==")
print(f"主终点（宏噬 DN↑, BH18）: submitter p.adj={mac_bh['p.adj_BH18']:.4g} 方向={mac_bh['direction']} | log2cpm p.adj={mac_bh_s['p.adj_BH18']:.4g} 方向={mac_bh_s['direction']}")
print(f"升格判据（主口径全样本 PDK2x宏噬 rho<0 & p<0.05）: submitter rho={mac_pri.rho} p={mac_pri.p} -> {'通过' if gate_primary else '未通过'}")
print(f"（敏感性口径同判据: log2cpm rho={mac_sens.rho} p={mac_sens.p} -> {'通过' if gate_sens else '未通过'}）")
print(f"最终判定: {'升格为第三免疫队列' if gate_primary else '混合/口径依赖 -> 必报规则：全表进补充表 + FFPE/试剂盒限定句，不升格 §3.7'}")
