# -*- coding: utf-8 -*-
"""
Module B1, second proteomics cohort: iProX IPX0003092000 / PXD040617
(glomerular iBAQ, DN n=53 vs control n=26; Zhao et al., J Transl Med 2023).

Zero handling (explicit): iBAQ values of 0 (deposited as 0.0 float; no empty
cells or NaN exist in this matrix) are treated as LEFT-CENSORED = not detected.
Zeros are NEVER imputed into medians. For every panel gene we report:
  - detection counts per group + two-sided Fisher exact test (detected/not);
  - among detected values only: median log2 iBAQ per group, log2FC = medDN-medC,
    two-sided Mann-Whitney U (computed only when >=3 detected values in BOTH groups);
  - a detection-differential direction (which group detects the gene more often).
Genes with <30% detection among controls are flagged detection-dominated
(effect estimates secondary; Fisher test primary). Genes absent from the
deposited matrix are included as rows with detected=0 and flag 'absent'.

Outputs (work folder):
  TableB1_iprox_panel.csv           deliverable panel table
  TableB1_iprox_global_context.csv  log2FC distribution context (>=50% detection both groups)
  iprox_analysis_summary.txt        console capture for report writing
"""
import os
import sys, io, re, json
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
import pandas as pd
from scipy import stats

BASE = BASE
CSV = BASE + r'\downloads_iprox\JTOM-2023-DN-Profiling.csv'

# ---------------------------------------------------------------- load
df = pd.read_csv(CSV)
GENE = 'GeneSymbol'
dn_cols = [c for c in df.columns if re.match(r'^DN\d+_iBAQ$', c)]
c_cols = [c for c in df.columns if re.match(r'^C\d+_iBAQ$', c)]
assert len(dn_cols) == 53 and len(c_cols) == 26, (len(dn_cols), len(c_cols))
N_DN, N_C = len(dn_cols), len(c_cols)
print(f'matrix: {df.shape[0]} genes x {N_DN} DN + {N_C} C iBAQ columns; '
      f'0-cells {(df[dn_cols+c_cols]==0).sum().sum()}, NaN {df[dn_cols+c_cols].isna().sum().sum()}')

V = df.set_index(GENE)[dn_cols + c_cols].astype(float)
LOG2 = np.log2(V.replace(0, np.nan))          # log2 of detected values only; 0 -> NaN (censored)

# ---------------------------------------------------------------- helpers
def gene_stats(gene, mat=LOG2):
    """Return stats dict for one gene symbol against matrix mat (log2, NaN=censored)."""
    out = {'gene': gene, 'detDN_n': f'0/{N_DN}', 'detC_n': f'0/{N_C}',
           'fisher_p': '', 'detection_direction': '', 'medDN': '', 'medC': '',
           'log2FC_detected_only': '', 'MWU_p': '', 'pct_log2FC': '',
           'interpretation_flag': 'absent'}
    if gene not in mat.index:
        return out
    x = mat.loc[gene, dn_cols].dropna()   # detected log2 values, DN
    y = mat.loc[gene, c_cols].dropna()    # detected log2 values, C
    detDN, detC = len(x), len(y)
    out['detDN_n'] = f'{detDN}/{N_DN}'
    out['detC_n'] = f'{detC}/{N_C}'
    if detDN == 0 and detC == 0:
        out['interpretation_flag'] = 'not_detected'
        return out
    # Fisher exact on detection table (two-sided)
    table = [[detDN, N_DN - detDN], [detC, N_C - detC]]
    orv, fp = stats.fisher_exact(table, alternative='two-sided')
    out['fisher_p'] = fp
    rateDN, rateC = detDN / N_DN, detC / N_C
    out['detection_direction'] = ('more_in_DN' if rateDN > rateC else
                                  'more_in_C' if rateC > rateDN else 'equal')
    # detected-only abundance stats
    if detDN:
        out['medDN'] = x.median()
    if detC:
        out['medC'] = y.median()
    if detDN and detC:
        out['log2FC_detected_only'] = x.median() - y.median()
    if detDN >= 3 and detC >= 3:
        out['MWU_p'] = stats.mannwhitneyu(x, y, alternative='two-sided').pvalue
    # flags
    flag = []
    if rateC < 0.30:
        flag.append('detection-dominated (sparse control, <30% of controls detected)')
    elif rateDN < 0.30:
        flag.append('detection-dominated (sparse DN, <30% of DN detected)')
    else:
        flag.append('quantified')
    if not (detDN >= 3 and detC >= 3):
        flag.append('MWU not computed (<3 detected values in a group)')
    elif min(detDN, detC) < 8:
        flag.append(f'low-n MWU ({min(detDN, detC)} vs {max(detDN, detC)})')
    out['interpretation_flag'] = '; '.join(flag)
    return out

# ---------------------------------------------------------------- global context
detDN_all = (V[dn_cols] > 0).sum(axis=1)
detC_all = (V[c_cols] > 0).sum(axis=1)
elig_mask = (detDN_all >= N_DN/2) & (detC_all >= N_C/2)
ggenes = LOG2.index[elig_mask]
gl2fc = pd.Series({g: LOG2.loc[g, dn_cols].median() - LOG2.loc[g, c_cols].median()
                   for g in ggenes}, name='log2FC')
n_dec = int((gl2fc < 0).sum()); n_inc = int((gl2fc > 0).sum())
print(f'\nglobal context: {len(gl2fc)} genes with >=50% detection in BOTH groups '
      f'(detDN>={int(np.ceil(N_DN/2))}, detC>={int(np.ceil(N_C/2))}); '
      f'decreased {n_dec} ({100*n_dec/len(gl2fc):.1f}%), increased {n_inc}; '
      f'median log2FC {gl2fc.median():+.3f} (IQR {gl2fc.quantile(.25):+.3f} to {gl2fc.quantile(.75):+.3f})')

pct = lambda v: stats.percentileofscore(gl2fc, v, kind='mean') if pd.notna(v) else np.nan

gctx = pd.DataFrame({
    'gene': gl2fc.index,
    'detDN_n': [f'{detDN_all[g]}/{N_DN}' for g in gl2fc.index],
    'detC_n': [f'{detC_all[g]}/{N_C}' for g in gl2fc.index],
    'log2FC_detected_only': gl2fc.values,
    'percentile_log2FC': [stats.percentileofscore(gl2fc, v, kind='mean') for v in gl2fc.values],
    'MWU_p': [stats.mannwhitneyu(LOG2.loc[g, dn_cols].dropna(), LOG2.loc[g, c_cols].dropna(),
                                 alternative='two-sided').pvalue for g in gl2fc.index],
})
gctx.to_csv(BASE + r'\TableB1_iprox_global_context.csv', index=False)

# ---------------------------------------------------------------- panel definition
# 36-node lactate-axis panel from TableB1_lactate_axis_panel.csv + module candidates
# (PTGS2, NFE2L2, MMP7) + face-validity / podocyte anchors.
axis_panel = ['LDHA','LDHB','SLC16A1','SLC16A3','SLC5A12','PDK1','PDK2','PDK3','PDK4',
              'PDHA1','PDHB','DLAT','DLD','MPC1','MPC2','PFKFB3','PFKL','PFKM','PFKP',
              'PKM','PKLR','ALDOB','ACLY','ACSS2','EP300','HDAC1','HDAC2','HDAC3',
              'SIRT1','SIRT2','SIRT3','SIRT5','ME1','ME2','KRT18','HIST1H4B']
candidates = ['PTGS2', 'NFE2L2', 'MMP7']
anchors = ['POSTN', 'VIM', 'COL1A1', 'COL3A1', 'NPHS1', 'NPHS2', 'PODXL', 'SYNPO']
panel = axis_panel + candidates + anchors

rows = []
for gene in panel:
    s = gene_stats(gene)
    # add global percentile for genes with an effect estimate (where it would sit)
    if s['log2FC_detected_only'] != '':
        s['pct_log2FC'] = pct(s['log2FC_detected_only'])
    rows.append(s)
res = pd.DataFrame(rows)

# order columns; numeric coercion for tidy output
num_cols = ['fisher_p', 'medDN', 'medC', 'log2FC_detected_only', 'MWU_p', 'pct_log2FC']
for c_ in num_cols:
    res[c_] = pd.to_numeric(res[c_], errors='coerce')
res = res[['gene', 'detDN_n', 'detC_n', 'fisher_p', 'detection_direction', 'medDN', 'medC',
           'log2FC_detected_only', 'MWU_p', 'pct_log2FC', 'interpretation_flag']]
res.to_csv(BASE + r'\TableB1_iprox_panel.csv', index=False, float_format='%.6g')

# ---------------------------------------------------------------- focused printout
print('\n=== PANEL TABLE (all rows) ===')
print(res.to_string(index=False))

print('\n=== PDH COMPLEX REPLICATION TEST (vs PXD041884 FFPE) ===')
ffpe = {'PDHA1': (-0.421486, 0.202020), 'PDHB': (-0.371667, 0.0176768),
        'DLAT': (-0.738347, 0.00505051), 'DLD': (-0.534047, 0.0176768)}
for g in ['PDHA1', 'PDHB', 'DLAT', 'DLD']:
    r = res[res.gene == g].iloc[0]
    f, fp = ffpe[g]
    print(f'{g}: iProX det {r.detDN_n} DN vs {r.detC_n} C | Fisher p={r.fisher_p:.4g} '
          f'| log2FC(det-only)={r.log2FC_detected_only:+.3f} | MWU p={r.MWU_p if pd.notna(r.MWU_p) else float("nan"):.4g} '
          f'| pct={r.pct_log2FC:.1f} || FFPE log2FC={f:+.3f} P={fp:.3g}')

print('\n=== MPC2 / MPC1 / transporters focus ===')
for g in ['MPC2', 'MPC1', 'SLC16A3', 'SLC16A1']:
    r = res[res.gene == g].iloc[0]
    print(f'{g}: {r.detDN_n} DN, {r.detC_n} C, Fisher p={r.fisher_p if pd.notna(r.fisher_p) else "NA"}, '
          f'dir={r.detection_direction}, log2FC={r.log2FC_detected_only if r.log2FC_detected_only != "" else "NA"}, '
          f'flag={r.interpretation_flag}')
# MPC2 detected-only medians detail
for g in ['MPC2']:
    if g in LOG2.index:
        x = LOG2.loc[g, dn_cols].dropna(); y = LOG2.loc[g, c_cols].dropna()
        print(f'MPC2 detected values: DN med={x.median():.3f} (n={len(x)}), range {x.min():.2f}-{x.max():.2f}; '
              f'C med={y.median():.3f} (n={len(y)}), range {y.min():.2f}-{y.max():.2f}')

print('\n=== KRT18 (curated lactylation-machinery candidate) ===')
r = res[res.gene == 'KRT18'].iloc[0]
print(r.to_string())

print('\n=== FACE VALIDITY ===')
for g in ['POSTN', 'VIM', 'COL1A1', 'COL3A1', 'NPHS1', 'NPHS2', 'PODXL', 'SYNPO']:
    r = res[res.gene == g].iloc[0]
    print(f'{g}: det {r.detDN_n} vs {r.detC_n} | log2FC={r.log2FC_detected_only} | '
          f'MWU p={r.MWU_p} | Fisher p={r.fisher_p} | pct={r.pct_log2FC}')

print('\n=== WHERE PDH/LDH SIT IN GLOBAL DISTRIBUTION (percentile of log2FC) ===')
for g in ['PDHA1','PDHB','DLAT','DLD','LDHA','LDHB','PKM','PKLR','ALDOB','MPC2','KRT18','ACLY','ME1','ME2','ACSS2']:
    r = res[res.gene == g].iloc[0]
    in_g = g in gl2fc.index
    print(f'{g}: log2FC={r.log2FC_detected_only} pct={r.pct_log2FC} in_global_set={in_g} MWU_p={r.MWU_p} Fisher_p={r.fisher_p}')

# ---------------------------------------------------------------- sensitivity:
# per-sample median-centering normalization (each sample shifted by median of its
# detected log2 values), then recompute panel medians/log2FC -> does direction hold?
med_shift = LOG2.median(axis=0)
LOG2N = LOG2.sub(med_shift, axis=1)
print('\n=== SENSITIVITY: per-sample median-centering (detected values) ===')
print('sample median shifts (log2): DN mean %+.3f sd %.3f; C mean %+.3f sd %.3f' % (
    med_shift[dn_cols].mean(), med_shift[dn_cols].std(), med_shift[c_cols].mean(), med_shift[c_cols].std()))
for g in ['PDHA1','PDHB','DLAT','DLD','MPC2','LDHA','LDHB','PKM','ALDOB','KRT18','POSTN','VIM']:
    if g in LOG2N.index:
        x = LOG2N.loc[g, dn_cols].dropna(); y = LOG2N.loc[g, c_cols].dropna()
        if len(x) and len(y):
            print(f'{g}: normalized log2FC={x.median()-y.median():+.3f} (raw {res[res.gene==g].iloc[0].log2FC_detected_only:+.3f})')

# ---------------------------------------------------------------- detection-background calibration
print('\n=== DETECTION BACKGROUND CALIBRATION ===')
det_per_sample_DN = (V[dn_cols] > 0).sum(axis=0)
det_per_sample_C = (V[c_cols] > 0).sum(axis=0)
print(f'genes detected per sample: DN mean {det_per_sample_DN.mean():.0f} '
      f'(sd {det_per_sample_DN.std():.0f}, range {det_per_sample_DN.min()}-{det_per_sample_DN.max()}), '
      f'C mean {det_per_sample_C.mean():.0f} (sd {det_per_sample_C.std():.0f}, '
      f'range {det_per_sample_C.min()}-{det_per_sample_C.max()}), MWU p='
      f'{stats.mannwhitneyu(det_per_sample_DN, det_per_sample_C, alternative="two-sided").pvalue:.3g}')
# genome-wide Fisher detection tests: how common is more_in_DN vs more_in_C at p<0.05?
gw_dir, gw_p = [], []
for g_ in V.index:
    a = int((V.loc[g_, dn_cols] > 0).sum()); b = int((V.loc[g_, c_cols] > 0).sum())
    if a + b == 0 or (a == 0 and b == 0):
        continue
    _, p_ = stats.fisher_exact([[a, N_DN - a], [b, N_C - b]], alternative='two-sided')
    gw_p.append(p_)
    r1, r2 = a / N_DN, b / N_C
    gw_dir.append('DN' if r1 > r2 else 'C' if r2 > r1 else 'eq')
import collections
cnt = collections.Counter(gw_dir)
sig = [(p_, d_) for p_, d_ in zip(gw_p, gw_dir) if p_ < 0.05]
sigc = collections.Counter(d_ for _, d_ in sig)
print(f'genes with any detection: {len(gw_p)}; direction more_in_DN {cnt["DN"]}, more_in_C {cnt["C"]}, equal {cnt["eq"]}')
print(f'Fisher detection p<0.05: {len(sig)} genes ({100*len(sig)/len(gw_p):.1f}%), of which more_in_DN {sigc["DN"]}, more_in_C {sigc["C"]}')
print(f'(expected at 5% alpha by chance: ~{0.05*len(gw_p):.0f})')
# histone H4 symbol check
h4 = [g_ for g_ in V.index if re.match(r'^H(IST1)?H4|^H4C', g_)]
print('histone H4-family symbols in matrix:', h4)
for g_ in h4:
    s_ = gene_stats(g_)
    print(f'  {g_}: {s_["detDN_n"]} vs {s_["detC_n"]}, log2FC={s_["log2FC_detected_only"]}, MWU p={s_["MWU_p"]}')

# ---------------------------------------------------------------- manifest digest
xml = open(BASE + r'\downloads_iprox\iprox_IPX0003092000.xml', encoding='utf-8', errors='replace').read()
print('\n=== MANIFEST DIGEST (iprox_IPX0003092000.xml) ===')
for pat, lab in [(r'PXD\d+', 'PXD'), (r'10\.\d{4,5}/[^\s"<]+', 'DOI'),
                 (r'<[^>]*name[^>]*>([^<]{1,120})', 'name-tags'),
                 (r'(JTOM[^<\s]*)', 'file-tags'), (r'(iBAQ)', 'iBAQ-mentions'),
                 (r'(Q Exactive[^<\n]{0,40})', 'instrument'), (r'(\d{4}-\d{2}-\d{2})', 'dates')]:
    hits = re.findall(pat, xml)
    if hits:
        print(f'{lab}: {sorted(set(hits))[:8]}')
print('\nDONE.')
