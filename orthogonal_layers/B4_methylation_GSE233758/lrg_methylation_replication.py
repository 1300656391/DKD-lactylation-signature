"""B4 analysis: LRG-list sweep on replication GSE145745 (EPICv2 PTEC, untreated T2D n=4 vs ctrl n=4).

Same gene/region assignment as the primary sweep; probes restricted to those present
on the EPIC v2 matrix. Direction-concordance report vs the primary GSE233758 sweep.
Output: TableB4_lrg_replication.csv
"""
import os
import gzip, time
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

BASE = BASE
SM = BASE + r'\downloads\GSE145745_series_matrix.txt.gz'

t0 = time.time()
samples, descs = [], []
rows, colnames = [], []
header_done = False
with gzip.open(SM, 'rt', encoding='utf-8', errors='replace') as f:
    for line in f:
        if not header_done:
            if line.startswith('!Sample_geo_accession'):
                samples = [x.strip('"') for x in line.rstrip('\n').split('\t')[1:]]
            elif line.startswith('!Sample_description'):
                descs = [x.strip('"') for x in line.rstrip('\n').split('\t')[1:]]
            elif line.startswith('!series_matrix_table_begin'):
                header_done = True
                continue
        else:
            if line.startswith('!series_matrix_table_end'):
                break
            parts = line.rstrip('\n').split('\t')
            if parts[0].strip('"') == 'ID_REF':
                colnames = [x.strip('"') for x in parts[1:]]
                continue
            rows.append(parts)
print('parsed %.1f s; rows %d cols %d' % (time.time() - t0, len(rows), len(colnames)))

ss = pd.DataFrame({'gsm': samples, 'desc': descs})
unt = ss[ss.desc.str.endswith('_Control')]
gsm_t2d = unt.loc[unt.desc == 'T2D_Control', 'gsm'].tolist()
gsm_norm = unt.loc[unt.desc == 'Normal_Control', 'gsm'].tolist()
print('untreated T2D n=%d  Normal n=%d' % (len(gsm_t2d), len(gsm_norm)))

# primary-sweep probe->gene/region assignment (from TableB4_lrg_list_methylation source = anno parse)
beta_ann = pd.read_pickle(BASE + r'\anno_epic_full.pkl')
prim = pd.read_csv(BASE + r'\TableB4_lrg_list_methylation.csv')

PROMOTER = {"TSS1500", "TSS200", "5'UTR", "1stExon"}
BODY = {"Body", "3'UTR", "ExonBnd"}
RANK = {"TSS200": 0, "TSS1500": 1, "5'UTR": 2, "1stExon": 3, "Body": 4, "3'UTR": 5, "ExonBnd": 6}

import pandas as pd2
esm = pd.read_csv(os.path.join(BASE, 'ESM_4.csv'))
gene_map = {}
for _, r in esm.iterrows():
    gene_map[str(r.symbol).strip()] = r.symbol
    al = str(r.alias).strip() if pd.notna(r.alias) else ''
    if al:
        gene_map[al] = r.symbol
LEGACY = {'AARS1': 'AARS', 'LARS1': 'LARS', 'MPC1': 'MPC1'}
names_all = set(gene_map) | set(LEGACY.values())

hits = []
for pid, nm, gr in zip(beta_ann.index, beta_ann['UCSC_RefGene_Name'].fillna(''),
                       beta_ann['UCSC_RefGene_Group'].fillna('')):
    if nm == '':
        continue
    ns = str(nm).split(';')
    gs = str(gr).split(';')
    if len(gs) == 1 and len(ns) > 1:
        gs = gs * len(ns)
    for n, g in zip(ns, gs):
        if n in names_all:
            sym = gene_map.get(n, n)
            cat = 'promoter' if g in PROMOTER else ('body' if g in BODY else None)
            if cat:
                hits.append((pid, sym, cat, RANK.get(g, 99)))
ha = pd.DataFrame(hits, columns=['probe', 'gene', 'cat', 'rank']).sort_values('rank') \
        .drop_duplicates(['probe', 'gene'], keep='first')
candset = set(ha.probe)
print('LRG probes in v1 annotation:', len(candset))

ids = [r[0].strip('"') for r in rows]
keep = [i for i, pid in enumerate(ids) if pid in candset]
print('LRG probes present on EPICv2 matrix:', len(keep))

mat = np.full((len(keep), len(colnames)), np.nan)
for ri, row_i in enumerate(keep):
    r = rows[row_i]
    for j in range(len(colnames)):
        v = r[j + 1].strip('"')
        if v.upper() not in ('NA', '', 'NULL'):
            mat[ri, j] = float(v)
beta = pd.DataFrame(mat, index=[ids[i] for i in keep], columns=samples)
print('NA cells:', int(np.isnan(mat).sum()), '/', mat.size)

res = []
for (gene, cat), g in ha[ha.probe.isin(beta.index)].groupby(['gene', 'cat']):
    probes = g.probe.tolist()
    sub = beta.loc[probes]
    dn = sub[gsm_t2d].mean(axis=0, skipna=True).dropna().values
    ct = sub[gsm_norm].mean(axis=0, skipna=True).dropna().values
    if len(dn) < 3 or len(ct) < 3:
        continue
    db = float(dn.mean() - ct.mean())
    p = float(mannwhitneyu(dn, ct, alternative='two-sided').pvalue)
    res.append(dict(gene=gene, region=cat, n_probes=len(probes), delta_beta_repl=db, wilcox_p_repl=p))
tab = pd.DataFrame(res)
for cat in ['promoter', 'body']:
    m = tab.region == cat
    ps = tab.loc[m, 'wilcox_p_repl'].values
    k = len(ps)
    o = np.argsort(ps)
    q = np.minimum.accumulate((ps[o] * k / np.arange(1, k + 1))[::-1])[::-1]
    bh = np.empty(k)
    bh[o] = np.clip(q, 0, 1)
    tab.loc[m, 'p_BH_repl'] = bh
tab.to_csv(BASE + r'\TableB4_lrg_replication.csv', index=False, float_format='%.6g')
print('WROTE TableB4_lrg_replication.csv rows:', len(tab))

# ---- direction concordance vs primary ----
mm = prim.merge(tab, on=['gene', 'region'], how='inner', suffixes=('', '_repl'))
both = mm.dropna(subset=['delta_beta', 'delta_beta_repl'])
print()
print('genes testable in BOTH cohorts (per region):')
for cat in ['promoter', 'body']:
    t = both[both.region == cat]
    same = (np.sign(t.delta_beta) == np.sign(t.delta_beta_repl))
    from scipy.stats import binomtest
    print(f'  {cat}: n={len(t)} | direction concordant {same.sum()}/{len(t)} '
          f'({100 * same.mean():.0f}%) | binom p={binomtest(int(same.sum()), len(t), 0.5).pvalue:.4f}')

print()
print('=== primary FDR<0.05 hits in replication ===')
for _, r in prim[prim.p_BH < 0.05].iterrows():
    q = tab[(tab.gene == r.gene) & (tab.region == r.region)]
    if len(q):
        q = q.iloc[0]
        agree = np.sign(r.delta_beta) == np.sign(q.delta_beta_repl)
        print(f'  {r.gene} {r.region}: primary db={r.delta_beta:+.4f} (FDR {r.p_BH:.3f}) | '
              f'repl db={q.delta_beta_repl:+.4f} p={q.wilcox_p_repl:.3f} | direction {"CONCORDANT" if agree else "DISCORDANT"}')
    else:
        print(f'  {r.gene} {r.region}: not testable in replication')

print()
print('=== PDK2/PDK4 in replication ===')
for gene in ['PDK2', 'PDK4']:
    for cat in ['promoter', 'body']:
        q = tab[(tab.gene == gene) & (tab.region == cat)]
        if len(q):
            q = q.iloc[0]
            pr = prim[(prim.gene == gene) & (prim.region == cat)].iloc[0]
            print(f'  {gene} {cat}: primary db={pr.delta_beta:+.4f} | repl db={q.delta_beta_repl:+.4f} '
                  f'p={q.wilcox_p_repl:.3f} | {"CONCORDANT" if np.sign(pr.delta_beta)==np.sign(q.delta_beta_repl) else "DISCORDANT"}')
