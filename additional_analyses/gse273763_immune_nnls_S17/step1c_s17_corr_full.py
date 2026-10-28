# evidence_gain/4_new_cohort/step1c_s17_corr_full.py ----
# P1A 建议采纳：S17 免疫必报表的 PDK x 18 细胞相关全表（双口径 x 全样本/DN 内，含 p 值）
import gzip, os
import numpy as np, pandas as pd
from scipy.stats import spearmanr

BASE = os.path.dirname(os.path.abspath(__file__))
CELLS = ["B_cells","T_cells","T_helper_cells","Tcm","Tem","Th1_cells","Th2_cells","TFH",
         "CD8_T_cells","Cytotoxic_cells","NK_cells","NK_CD56dim_cells","DC","iDC","Eosinophils",
         "Macrophages","Mast_cells","Neutrophils"]

rows = []
for tag in ("submitter", "log2cpm"):
    ss = pd.read_csv(os.path.join(BASE, f"gse273763_ssgsea_{tag}.csv")).set_index("sample")
    if tag == "submitter":
        with gzip.open(os.path.join(BASE, "GSE273763_DKD_GRCH38_forGEO_July2026_v2.txt.gz"), "rt", encoding="utf-8", errors="replace") as f:
            M = pd.read_csv(f, sep="\t", index_col=0)
    else:
        with gzip.open(os.path.join(BASE, "GSE273763_counts_symbol.csv.gz"), "rt", encoding="utf-8", errors="replace") as f:
            C = pd.read_csv(f, index_col=0)
        M = np.log2(C.div(C.sum(axis=0), axis=1) * 1e6 + 1)
    M.columns = M.columns.str.replace(r"\s+", "", regex=True)
    M.index = M.index.astype(str).str.upper()
    M = M.loc[~M.index.duplicated(keep="first")]
    ss = ss.loc[M.columns]
    is_dn = (ss["group"] == "DKD").values
    for gene in ("PDK2", "PDK4"):
        v = M.loc[gene].values
        for cell in CELLS:
            s = ss[cell].values
            for scope, mask in (("all", np.ones(len(ss), bool)), ("DN_only", is_dn)):
                r, p = spearmanr(v[mask], s[mask])
                rows.append(dict(caliber=tag, gene=gene, cell=cell, scope=scope,
                                 rho=round(r, 3), p=round(p, 6)))
out = pd.DataFrame(rows)
out.to_csv(os.path.join(BASE, "gse273763_pdk_corr_full.csv"), index=False)
print(f"rows: {len(out)} -> gse273763_pdk_corr_full.csv")
print(out[(out.cell.isin(['Macrophages', 'B_cells'])) & (out.gene == 'PDK2')].to_string(index=False))
print("S17_CORR_FULL_DONE")
