"""Build final TableB4_pdk_methylation.csv (per-probe rows + category summary rows)
plus a sensitivity pass excluding SNP-affected probes."""
import os
import numpy as np
import pandas as pd
from scipy import stats

BASE = BASE
tab = pd.read_csv(BASE + r'\probe_stats_raw.csv')
snps = pd.read_csv(BASE + r'\snp_flags.csv')
snpmap = dict(zip(snps['probe'], snps['snp_hits']))
tab['snp_flag'] = tab['cgID'].map(lambda p: snpmap.get(p, ''))

# ---- sensitivity: exclude SNP-affected probes, re-BH ----
clean = tab[tab['snp_flag'] == ''].copy()
k2 = len(clean)
o = np.argsort(clean['p'].values)
sp = clean['p'].values[o]
q = np.minimum.accumulate((sp * k2 / np.arange(1, k2 + 1))[::-1])[::-1]
bh2 = np.empty(k2)
bh2[o] = np.clip(q, 0, 1)
clean['p_BH'] = bh2
print('sensitivity family k =', k2)
for (gene, cat), g in clean.groupby(['gene', 'region_category']):
    print(f'[sensitivity, SNP-clean] {gene} {cat}: k={len(g)} mean_db={g.delta_beta.mean():+.4f} '
          f'range=[{g.delta_beta.min():+.4f},{g.delta_beta.max():+.4f}] min_p={g.p.min():.4g} min_BH={g.p_BH.min():.4g}')
clean.to_csv(BASE + r'\probe_stats_snp_clean.csv', index=False)

# ---- per-probe table (required columns first, extras after) ----
cols = ['cgID', 'gene', 'region_category', 'mean_beta_DN', 'mean_beta_ctrl', 'delta_beta',
        'Mvalue_logFC', 'p', 'p_BH', 'refgene_group', 'chr', 'pos', 'relation_to_island',
        'snp_flag', 'p_batch_adj', 'n_DN', 'n_ctrl']
out = tab[cols].copy().sort_values(['gene', 'region_category', 'p'])

def fmt(v, nd=4):
    return round(float(v), nd)

summary_rows = []
for (gene, cat), g in tab.groupby(['gene', 'region_category']):
    summary_rows.append({
        'cgID': f'SUMMARY {gene} {cat}',
        'gene': gene, 'region_category': cat,
        'mean_beta_DN': g.mean_beta_DN.mean(), 'mean_beta_ctrl': g.mean_beta_ctrl.mean(),
        'delta_beta': g.delta_beta.mean(),
        'Mvalue_logFC': g.Mvalue_logFC.mean(),
        'p': g.p.min(), 'p_BH': g.p_BH.min(),
        'refgene_group': f'k={len(g)}; range delta_beta=[{g.delta_beta.min():+.4f},{g.delta_beta.max():+.4f}]; '
                         f'n_probes_delta_beta>0={(g.delta_beta > 0).sum()}/{len(g)}',
        'chr': '', 'pos': '', 'relation_to_island': '',
        'snp_flag': f'SNP-clean k={(g.snp_flag == "").sum()}',
        'p_batch_adj': g.p_batch_adj.min(),
        'n_DN': '', 'n_ctrl': '',
    })
summ = pd.DataFrame(summary_rows)
final = pd.concat([out, summ], ignore_index=True)
final.to_csv(BASE + r'\TableB4_pdk_methylation.csv', index=False)
print('wrote TableB4_pdk_methylation.csv rows:', len(final))

# ---- region-level tests to a small stats file for the report ----
beta = pd.read_pickle(BASE + r'\beta_matrix.pkl')
ss = pd.read_csv(BASE + r'\sample_sheet.csv')
gsm_dn = ss.loc[ss.group == 'DN', 'gsm'].tolist()
gsm_ct = ss.loc[ss.group == 'Control', 'gsm'].tolist()
batch = ss.set_index('gsm')['batch']
EPS = 1e-3
def mval(b):
    return np.log2((b + EPS) / (1 - b + EPS))

region_stats = []
def ols_ttest(y, X, col):
    XtX_inv = np.linalg.inv(X.T @ X)
    beta_hat = XtX_inv @ X.T @ y
    resid = y - X @ beta_hat
    dfres = X.shape[0] - X.shape[1]
    s2 = resid @ resid / dfres
    se = np.sqrt(s2 * XtX_inv[col, col])
    tval = beta_hat[col] / se
    p = 2 * stats.t.sf(abs(tval), dfres)
    return beta_hat[col], tval, p

for gene in ('PDK2', 'PDK4'):
    for cat in ('promoter', 'body'):
        for label, sub in [('all', tab), ('snp_clean', clean)]:
            probes = sub.loc[(sub.gene == gene) & (sub.region_category == cat), 'cgID']
            if len(probes) == 0:
                continue
            mm = mval(beta.loc[probes].values)
            reg = pd.Series(np.nanmean(mm, axis=0), index=beta.columns)
            # unadjusted
            dn = reg[gsm_dn].dropna().values; ct = reg[gsm_ct].dropna().values
            tt = stats.ttest_ind(dn, ct)
            # batch adjusted
            keep = [g for g in gsm_ct + gsm_dn if not pd.isna(reg[g])]
            y2 = np.array([reg[g] for g in keep])
            X2 = np.column_stack([np.ones(len(keep)),
                                  np.array([1.0 if g in gsm_dn else 0.0 for g in keep]),
                                  np.array([1.0 if batch[g] == 'B2' else 0.0 for g in keep])])
            c2, t2, p2 = ols_ttest(y2, X2, 1)
            db = beta.loc[probes, gsm_dn].mean(axis=0).mean() - beta.loc[probes, gsm_ct].mean(axis=0).mean()
            region_stats.append(dict(gene=gene, region=cat, set=label, k=len(probes),
                                     mean_delta_beta=db, t=tt.statistic, p=tt.pvalue,
                                     t_batch_adj=t2, p_batch_adj=p2))
rs = pd.DataFrame(region_stats)
rs.to_csv(BASE + r'\region_level_stats.csv', index=False)
print(rs.to_string(index=False))
