"""B4: PDK2/PDK4 promoter methylation analysis in GSE233758 (DN vs control PT)."""
import os
import numpy as np
import pandas as pd
from scipy import stats

BASE = BASE

beta = pd.read_pickle(BASE + r'\beta_matrix.pkl')
ss = pd.read_csv(BASE + r'\sample_sheet.csv')
anno = pd.read_pickle(BASE + r'\anno_epic_full.pkl')

# ---------- probe set: PDK2 / PDK4 ----------
PROMOTER = {"TSS1500", "TSS200", "5'UTR", "1stExon"}
BODY = {"Body", "3'UTR", "ExonBnd"}
RANK = {"TSS200": 0, "TSS1500": 1, "5'UTR": 2, "1stExon": 3, "Body": 4, "3'UTR": 5, "ExonBnd": 6}

rows = []
mask = anno['UCSC_RefGene_Name'].fillna('').str.contains('PDK2|PDK4', regex=True)
for pid, r in anno[mask].iterrows():
    names = str(r['UCSC_RefGene_Name']).split(';')
    groups = str(r['UCSC_RefGene_Group']).split(';')
    if len(groups) == 1 and len(names) > 1:
        groups = groups * len(names)
    for gene in ('PDK2', 'PDK4'):
        hits = [g for n, g in zip(names, groups) if n == gene]
        if hits:
            best = sorted(set(hits), key=lambda g: RANK.get(g, 99))[0]
            cat = 'promoter' if best in PROMOTER else ('body' if best in BODY else 'other')
            rows.append(dict(probe=pid, gene=gene, best_group=best, region_category=cat,
                             all_groups=';'.join(sorted(set(hits))),
                             chr=r['chr'], pos=int(r['pos']),
                             relation_to_island=r['Relation_to_Island']))
cand = pd.DataFrame(rows).set_index('probe')
print('candidate probes:', len(cand), 'by gene/category:')
print(cand.groupby(['gene', 'region_category']).size())

# present in beta matrix?
cand = cand.loc[cand.index.intersection(beta.index)]
print('present in matrix:', len(cand))

# ---------- groups ----------
gsm_dn = ss.loc[ss.group == 'DN', 'gsm'].tolist()
gsm_ct = ss.loc[ss.group == 'Control', 'gsm'].tolist()
batch = ss.set_index('gsm')['batch']
print('n DN =', len(gsm_dn), ' n Ctrl =', len(gsm_ct))

EPS = 1e-3

def mval(b):
    return np.log2((b + EPS) / (1 - b + EPS))

def ols_ttest(y, X, col):
    """OLS; returns coef, se, t, p for column col."""
    XtX_inv = np.linalg.inv(X.T @ X)
    beta_hat = XtX_inv @ X.T @ y
    resid = y - X @ beta_hat
    dfres = X.shape[0] - X.shape[1]
    s2 = resid @ resid / dfres
    se = np.sqrt(s2 * XtX_inv[col, col])
    tval = beta_hat[col] / se
    p = 2 * stats.t.sf(abs(tval), dfres)
    return beta_hat[col], se, tval, p

res = []
for pid, r in cand.iterrows():
    b = beta.loc[pid]
    bd = b[gsm_dn]
    bc = b[gsm_ct]
    nd, nc = bd.notna().sum(), bc.notna().sum()
    if nd < 10 or nc < 10:
        print('skip (too few values)', pid, nd, nc)
        continue
    m = mval(b.values)
    m = pd.Series(m, index=b.index)
    dn = m[gsm_dn].dropna().values
    ct = m[gsm_ct].dropna().values
    # design: intercept + group(DN=1); then sensitivity: + batch(B2)
    y = np.concatenate([ct, dn])
    x1 = np.column_stack([np.ones(len(y)), np.r_[np.zeros(len(ct)), np.ones(len(dn))]])
    coef, se, tval, p = ols_ttest(y, x1, 1)
    # batch-adjusted (batch from sample sheet; only samples with non-NA M)
    keep = [g for g in gsm_ct + gsm_dn if not pd.isna(b[g])]
    yy = np.array([m[g] for g in keep])
    X2 = np.column_stack([
        np.ones(len(keep)),
        np.array([1.0 if g in gsm_dn else 0.0 for g in keep]),
        np.array([1.0 if batch[g] == 'B2' else 0.0 for g in keep]),
    ])
    coef2, se2, t2, p2 = ols_ttest(yy, X2, 1)
    res.append(dict(
        cgID=pid, gene=r['gene'], region_category=r['region_category'],
        refgene_group=r['best_group'], chr=r['chr'], pos=r['pos'],
        relation_to_island=r['relation_to_island'],
        mean_beta_DN=float(bd.mean()), mean_beta_ctrl=float(bc.mean()),
        delta_beta=float(bd.mean() - bc.mean()),
        Mvalue_logFC=float(coef), t=float(tval), p=float(p),
        p_batch_adj=float(p2), delta_beta_batch_adj=float(coef2),
        n_DN=int(nd), n_ctrl=int(nc),
    ))

tab = pd.DataFrame(res)
k = len(tab)
# BH across candidate family
o = np.argsort(tab['p'].values)
m = k
sorted_p = tab['p'].values[o]
q = np.minimum.accumulate((sorted_p * m / np.arange(1, m + 1))[::-1])[::-1]
bh = np.empty(m)
bh[o] = np.clip(q, 0, 1)
tab['p_BH'] = bh
print('family size k =', k)
tab.to_csv(BASE + r'\probe_stats_raw.csv', index=False)

pd.set_option('display.width', 250)
print(tab[['cgID', 'gene', 'region_category', 'refgene_group', 'mean_beta_DN', 'mean_beta_ctrl',
           'delta_beta', 'Mvalue_logFC', 'p', 'p_BH', 'p_batch_adj', 'n_DN', 'n_ctrl']].to_string(index=False))

# ---------- category summaries ----------
print('\n=== category summaries (delta_beta = DN - Ctrl; positive = hypermethylation in DN) ===')
for (gene, cat), g in tab.groupby(['gene', 'region_category']):
    print(f'{gene} {cat}: k={len(g)} mean_db={g.delta_beta.mean():+.4f} '
          f'range=[{g.delta_beta.min():+.4f}, {g.delta_beta.max():+.4f}] '
          f'min_p={g.p.min():.4g} min_BH={g.p_BH.min():.4g}')

# region-level test: average M across promoter probes per sample, then DN vs Ctrl
for gene in ('PDK2', 'PDK4'):
    for cat in ('promoter', 'body'):
        probes = tab.loc[(tab.gene == gene) & (tab.region_category == cat), 'cgID']
        if len(probes) == 0:
            continue
        mm = mval(beta.loc[probes].values)  # probes x samples
        reg = pd.Series(np.nanmean(mm, axis=0), index=beta.columns)
        dn = reg[gsm_dn].dropna().values
        ct = reg[gsm_ct].dropna().values
        tt = stats.ttest_ind(dn, ct)
        db = beta.loc[probes, gsm_dn].mean(axis=0).mean() - beta.loc[probes, gsm_ct].mean(axis=0).mean()
        print(f'{gene} {cat} region-level: k={len(probes)} mean_delta_beta={db:+.4f} '
              f't_stat={tt.statistic:+.3f} p={tt.pvalue:.4g}')

# overall direction sign test across promoter probes
for gene in ('PDK2', 'PDK4'):
    g = tab[(tab.gene == gene) & (tab.region_category == 'promoter')]
    npos = int((g.delta_beta > 0).sum())
    print(f'{gene} promoter: {npos}/{len(g)} probes positive delta-beta')
