"""B4 analysis: 171-gene lactylation-axis (LRG) list-wide methylation sweep on GSE233758.

Per gene x region class (promoter / gene body):
  - probes annotated to the gene (UCSC_RefGene via ilm10b4 hg19 anno; same parsing as s3)
  - per-sample mean beta across probes (>=1 probe)
  - delta-beta = mean(DN) - mean(ctrl); Wilcoxon rank-sum (exact, two-sided)
  - BH-FDR across genes WITHIN each region class
Output: TableB4_lrg_list_methylation.csv + summary printout.
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

BASE = BASE
ESM = os.path.join(BASE, 'ESM_4.csv')

beta = pd.read_pickle(BASE + r'\beta_matrix.pkl')
ss = pd.read_csv(BASE + r'\sample_sheet.csv')
anno = pd.read_pickle(BASE + r'\anno_epic_full.pkl')

gsm_dn = ss.loc[ss.group == 'DN', 'gsm'].tolist()
gsm_ct = ss.loc[ss.group == 'Control', 'gsm'].tolist()
print('n DN/ctrl:', len(gsm_dn), len(gsm_ct))

lrg = pd.read_csv(ESM)
gene_map = {}  # search-symbol -> ESM symbol
for _, r in lrg.iterrows():
    gene_map[str(r.symbol).strip()] = r.symbol
    al = str(r.alias).strip() if pd.notna(r.alias) else ''
    if al:
        gene_map[al] = r.symbol

# UCSC RefGene legacy-name fixes for symbols not present under current HGNC names
LEGACY = {
    'H3-3A': ['H3F3A'], 'H3C1': ['HIST1H3A'], 'H3C2': ['HIST1H3B'], 'H3C3': ['HIST1H3C'],
    'H4C1': ['HIST1H4A'], 'H4C2': ['HIST1H4B'], 'H4C11': ['HIST1H4K'], 'H4C12': ['HIST1H4L'],
    'H4C13': ['HIST1H4M'], 'H4C14': ['HIST1H4N'], 'H4C15': ['HIST1H4I'], 'H4C16': ['HIST1H4J'],
    'H2BC4': ['HIST1H2BC'], 'H2BC21': ['HIST2H2BE'], 'LDHAL6A': ['LDHAL6A'],
    'ADPGK': ['ADPGK'], 'MPC1': ['MPC1', 'BRP44L'], 'SLC5A12': ['SLC5A12'],
    'AARS1': ['AARS'], 'LARS1': ['LARS'],
}
for sym, legacy in LEGACY.items():
    if sym in set(lrg.symbol):
        for lg in legacy:
            gene_map[lg] = sym

names_all = set(gene_map)

PROMOTER = {"TSS1500", "TSS200", "5'UTR", "1stExon"}
BODY = {"Body", "3'UTR", "ExonBnd"}
RANK = {"TSS200": 0, "TSS1500": 1, "5'UTR": 2, "1stExon": 3, "Body": 4, "3'UTR": 5, "ExonBnd": 6}

# probe -> (esm_symbol, category) for every LRG gene hit (probe may serve 2 genes)
ann_names = anno['UCSC_RefGene_Name'].fillna('')
ann_groups = anno['UCSC_RefGene_Group'].fillna('')
hits = []  # (probe, esm_symbol, best_group)
for pid, nm, gr in zip(anno.index, ann_names, ann_groups):
    if nm == '':
        continue
    ns = str(nm).split(';')
    gs = str(gr).split(';')
    if len(gs) == 1 and len(ns) > 1:
        gs = gs * len(ns)
    for n, g in zip(ns, gs):
        if n in names_all:
            hits.append((pid, gene_map[n], g))
hit_df = pd.DataFrame(hits, columns=['probe', 'gene', 'group'])
print('probe-gene assignments:', len(hit_df), '| LRG genes with >=1 annotated probe:',
      hit_df.gene.nunique(), 'of', lrg.symbol.nunique())

# dedupe (probe, gene): keep promoter-est group
hit_df['rank'] = hit_df.group.map(RANK).fillna(99)
hit_df = (hit_df.sort_values('rank')
                .drop_duplicates(['probe', 'gene'], keep='first'))

def best_cat(g):
    return 'promoter' if g in PROMOTER else ('body' if g in BODY else 'other')

hit_df['category'] = hit_df.group.map(best_cat)
hit_df = hit_df[hit_df.category != 'other']

present = hit_df[hit_df.probe.isin(beta.index)]
print('assignments present in matrix:', len(present), 'genes:', present.gene.nunique())

missing = sorted(set(lrg.symbol) - set(present.gene))
print('LRG genes with NO annotated+present probe (%d):' % len(missing), ', '.join(missing))

rows = []
for (gene, cat), g in present.groupby(['gene', 'category']):
    probes = g.probe.tolist()
    sub = beta.loc[probes]                     # probes x 33
    samp_mean = sub.mean(axis=0, skipna=True)  # per-sample mean beta
    dn = samp_mean[gsm_dn].dropna().values
    ct = samp_mean[gsm_ct].dropna().values
    if len(dn) < 5 or len(ct) < 5:
        continue
    db = float(dn.mean() - ct.mean())
    p = float(mannwhitneyu(dn, ct, alternative='two-sided').pvalue)
    dn_m = sub[gsm_dn].mean(axis=1, skipna=True)
    ct_m = sub[gsm_ct].mean(axis=1, skipna=True)
    rows.append(dict(gene=gene, region=cat, n_probes=len(probes),
                     mean_beta_DN=float(dn.mean()), mean_beta_ctrl=float(ct.mean()),
                     delta_beta=db, wilcox_p=p,
                     n_pos_probes=int((dn_m - ct_m > 0).sum())))
tab = pd.DataFrame(rows)

# BH within region class
for cat, m in tab.groupby('region').groups.items():
    idx = m
    ps = tab.loc[idx, 'wilcox_p'].values
    k = len(ps)
    o = np.argsort(ps)
    q = np.minimum.accumulate((ps[o] * k / np.arange(1, k + 1))[::-1])[::-1]
    bh = np.empty(k)
    bh[o] = np.clip(q, 0, 1)
    tab.loc[idx, 'p_BH'] = bh

# rank of PDK2/PDK4
for cat in ['promoter', 'body']:
    t = tab[tab.region == cat].sort_values('wilcox_p').reset_index(drop=True)
    t['rank'] = t.index + 1
    for gene in ['PDK2', 'PDK4']:
        r = t[t.gene == gene]
        if len(r):
            print(f'rank {gene} {cat}: {int(r["rank"].iloc[0])}/{len(t)} (p={r.wilcox_p.iloc[0]:.4g}, '
                  f'db={r.delta_beta.iloc[0]:+.4f}, k={int(r.n_probes.iloc[0])})')
        else:
            print(f'{gene} {cat}: not in table')

tab['in_LRG_list'] = True
tab = tab.sort_values(['region', 'wilcox_p'])
tab.to_csv(BASE + r'\TableB4_lrg_list_methylation.csv', index=False, float_format='%.6g')
print('WROTE TableB4_lrg_list_methylation.csv rows:', len(tab))

pd.set_option('display.width', 240)
print()
for cat in ['promoter', 'body']:
    t = tab[tab.region == cat]
    sig = t[t.p_BH < 0.05]
    print(f'=== {cat}: {len(t)} genes testable | FDR<0.05: {len(sig)} '
          f'(hyper {int((sig.delta_beta > 0).sum())} / hypo {int((sig.delta_beta < 0).sum())})')
    print(f'    median delta-beta across axis genes: {t.delta_beta.median():+.4f} | '
          f'positive: {int((t.delta_beta > 0).sum())}/{len(t)} | p<0.05: {int((t.wilcox_p < 0.05).sum())}')
    print(sig.head(15)[['gene', 'n_probes', 'delta_beta', 'wilcox_p', 'p_BH', 'n_pos_probes']].to_string(index=False))
    print()

# metabolic-gene subset: does body hypermethylation concentrate in metabolic genes?
METABOLIC = set('''LDHA LDHB SLC16A1 SLC16A3 SLC5A8 SLC5A12 SLC2A1 SLC2A2 SLC2A3 PDK1 PDK2 PDK3 PDK4
PDHA1 PDHB PDHX DLAT DLD MPC1 MPC2 PCK1 PCK2 PC PKM PKLR PFKL PFKM PFKP PFKFB3 HK1 HK2
ALDOA ALDOB ALDOC GAPDH PGK1 ENO1 ENO2 ENO3 TPI1 CS ACO1 ACO2 IDH1 IDH2 IDH3A IDH3B IDH3G
OGDH DLST SUCLA2 SUCLG1 SUCLG2 SDHA SDHB SDHC SDHD FH MDH1 MDH2 ACLY ACSS1 ACSS2 ME1 ME2 ME3
GYS1 G6PC FBP1 PGM1 PGM2 ENO4 ADPGK HKDC1 LDHAL6A SLC16A7 SLC16A8 SLC5A6'''.split())
tab['metabolic'] = tab.gene.isin(METABOLIC)
for cat in ['promoter', 'body']:
    t = tab[tab.region == cat]
    a = t[t.metabolic].delta_beta
    b = t[~t.metabolic].delta_beta
    from scipy.stats import mannwhitneyu as mw
    pv = mw(a, b, alternative='two-sided').pvalue
    print(f'{cat}: metabolic-gene median db {a.median():+.4f} (n={len(a)}) vs other LRG {b.median():+.4f} '
          f'(n={len(b)}), MW p={pv:.3f}')
