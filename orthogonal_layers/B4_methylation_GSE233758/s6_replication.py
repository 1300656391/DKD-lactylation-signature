"""Replication: GSE145745 EPICv2 PTEC, untreated T2D vs untreated non-diabetic control."""
import os
import gzip
import time
import numpy as np
import pandas as pd
from scipy import stats

BASE = BASE
SM = BASE + r'\downloads\GSE145745_series_matrix.txt.gz'

def split_tabs(line):
    return line.rstrip('\n').split('\t')

titles, samples, descs = [], [], []
rows, colnames = [], []
header_done = False
t0 = time.time()
with gzip.open(SM, 'rt', encoding='utf-8', errors='replace') as f:
    for line in f:
        if not header_done:
            if line.startswith('!Sample_geo_accession'):
                samples = [x.strip('"') for x in split_tabs(line)[1:]]
            elif line.startswith('!Sample_title'):
                titles = [x.strip('"') for x in split_tabs(line)[1:]]
            elif line.startswith('!Sample_description'):
                descs = [x.strip('"') for x in split_tabs(line)[1:]]
            elif line.startswith('!series_matrix_table_begin'):
                header_done = True
                continue
        else:
            if line.startswith('!series_matrix_table_end'):
                break
            parts = split_tabs(line)
            if parts[0] == '"ID_REF"':
                colnames = [x.strip('"') for x in parts[1:]]
                continue
            rows.append(parts)
print('parsed', time.time() - t0, 's; rows:', len(rows), 'cols:', len(colnames))

# untreated samples only
ss = pd.DataFrame({'gsm': samples, 'title': titles, 'desc': descs})
unt = ss[ss.desc.str.endswith('_Control')]
gsm_t2d = unt.loc[unt.desc == 'T2D_Control', 'gsm'].tolist()
gsm_norm = unt.loc[unt.desc == 'Normal_Control', 'gsm'].tolist()
print('untreated T2D:', len(gsm_t2d), 'untreated Normal:', len(gsm_norm))

cand = pd.read_csv(BASE + r'\probe_stats_raw.csv')[['cgID', 'gene', 'region_category', 'refgene_group']]
candset = set(cand['cgID'])

ids = [r[0].strip('"') for r in rows]
keep_rows = [i for i, pid in enumerate(ids) if pid in candset]
print('candidate probes present in EPICv2 matrix:', len(keep_rows), 'of', len(candset))

mat = np.full((len(keep_rows), len(colnames)), np.nan)
for ri, row_i in enumerate(keep_rows):
    r = rows[row_i]
    for j in range(len(colnames)):
        v = r[j + 1].strip('"')
        if v not in ('NA', '', 'null'):
            mat[ri, j] = float(v)
beta = pd.DataFrame(mat, index=[ids[i] for i in keep_rows], columns=samples)
print('beta NA cells:', int(np.isnan(mat).sum()))

EPS = 1e-3
def mval(b):
    return np.log2((b + EPS) / (1 - b + EPS))

def ols_ttest(y, X, col):
    XtX_inv = np.linalg.inv(X.T @ X)
    beta_hat = XtX_inv @ X.T @ y
    resid = y - X @ beta_hat
    dfres = X.shape[0] - X.shape[1]
    s2 = resid @ resid / dfres
    se = np.sqrt(s2 * XtX_inv[col, col])
    tval = beta_hat[col] / se
    return beta_hat[col], tval, 2 * stats.t.sf(abs(tval), dfres)

res = []
for pid in beta.index:
    b = beta.loc[pid]
    bd, bc = b[gsm_t2d], b[gsm_norm]
    if bd.notna().sum() < 3 or bc.notna().sum() < 3:
        continue
    m = pd.Series(mval(b.values), index=b.index)
    dn = m[gsm_t2d].dropna().values
    ct = m[gsm_norm].dropna().values
    y = np.concatenate([ct, dn])
    X = np.column_stack([np.ones(len(y)), np.r_[np.zeros(len(ct)), np.ones(len(dn))]])
    coef, tval, p = ols_ttest(y, X, 1)
    meta = cand.set_index('cgID').loc[pid]
    res.append(dict(cgID=pid, gene=meta['gene'], region_category=meta['region_category'],
                    mean_beta_T2D=float(bd.mean()), mean_beta_ctrl=float(bc.mean()),
                    delta_beta=float(bd.mean() - bc.mean()), Mvalue_logFC=float(coef),
                    t=float(tval), p=float(p),
                    n_T2D=int(bd.notna().sum()), n_ctrl=int(bc.notna().sum())))

tab = pd.DataFrame(res)
k = len(tab)
o = np.argsort(tab['p'].values)
sp = tab['p'].values[o]
q = np.minimum.accumulate((sp * k / np.arange(1, k + 1))[::-1])[::-1]
bh = np.empty(k)
bh[o] = np.clip(q, 0, 1)
tab['p_BH'] = bh
tab.to_csv(BASE + r'\TableB4_replication.csv', index=False)
print('replication family k =', k)
pd.set_option('display.width', 220)
print(tab[['cgID', 'gene', 'region_category', 'mean_beta_T2D', 'mean_beta_ctrl', 'delta_beta',
           'Mvalue_logFC', 'p', 'p_BH']].to_string(index=False))

# region-level
for gene in ('PDK2', 'PDK4'):
    for cat in ('promoter', 'body'):
        probes = tab.loc[(tab.gene == gene) & (tab.region_category == cat), 'cgID']
        if len(probes) == 0:
            continue
        mm = mval(beta.loc[probes].values)
        reg = pd.Series(np.nanmean(mm, axis=0), index=beta.columns)
        dn = reg[gsm_t2d].dropna().values
        ct = reg[gsm_norm].dropna().values
        tt = stats.ttest_ind(dn, ct)
        db = beta.loc[probes, gsm_t2d].mean(axis=0).mean() - beta.loc[probes, gsm_norm].mean(axis=0).mean()
        npos = int((tab.loc[probes, 'delta_beta'] > 0).sum())
        print(f'{gene} {cat}: k={len(probes)} mean_db={db:+.4f} region_t={tt.statistic:+.2f} p={tt.pvalue:.4g} pos={npos}/{len(probes)}')
