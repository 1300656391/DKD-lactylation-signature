"""B1 module: PXD041884 FFPE kidney proteome, DKD vs control per-protein stats.

Input : downloads/Secocnd_Extraction_GPQ_Norm_TS_All_Peptides_FFPE_Protein_Rollup_Workbook.xlsx
        sheet 'All Protien Rolled Up' (normalized total-protein-rollup log2 quantities)
Groups: DKD  n=5 (Sample_01_D1, 05_D3, 06_D4, 09_D5, 10_D6)
        Ctrl n=7 (Sample_03_C1, 04_C2d, 07_C3, 08_C4, 12_C6, 13_C7, 14_C9)
Stats : log2FC = median(log2 DKD) - median(log2 Ctrl); exact Wilcoxon rank-sum
        (Mann-Whitney U, two-sided, exact where no ties).
"""
import os
import openpyxl, os, json
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

BASE = os.path.dirname(os.path.abspath(__file__))
DL = os.path.join(BASE, 'downloads')
WORK = BASE

DKD_SAMPLES = ['Sample_01_D1', 'Sample_05_D3', 'Sample_06_D4', 'Sample_09_D5', 'Sample_10_D6']
CTRL_SAMPLES = ['Sample_03_C1', 'Sample_04_C2d', 'Sample_07_C3', 'Sample_08_C4',
                'Sample_12_C6', 'Sample_13_C7', 'Sample_14_C9']

fn = os.path.join(DL, 'Secocnd_Extraction_GPQ_Norm_TS_All_Peptides_FFPE_Protein_Rollup_Workbook.xlsx')
wb = openpyxl.load_workbook(fn, read_only=True, data_only=True)
ws = wb['All Protien Rolled Up']
rows = ws.iter_rows(values_only=True)
hdr = list(next(rows))

# map header -> column index (log2 columns are '<sample>: Log2')
colmap = {}
for i, h in enumerate(hdr):
    if h is None:
        continue
    h = str(h).strip()
    if h.endswith(': Log2'):
        colmap[h[:-len(': Log2')]] = i
    elif h == 'Sample_01_D1':  # duplicate unsuffixed column -> skip (flag column, not quant)
        pass

missing = [s for s in DKD_SAMPLES + CTRL_SAMPLES if s not in colmap]
assert not missing, f'missing sample columns: {missing}'

recs = []
for ri, row in enumerate(rows, start=2):
    name = str(row[0]) if row[0] is not None else ''
    acc_field = str(row[1]).strip() if row[1] is not None else ''
    peptide_count = row[4]
    primary = acc_field.split('|')[0]
    members = [a for a in acc_field.split('|') if a]
    dkd = np.array([row[colmap[s]] for s in DKD_SAMPLES], dtype=float)
    ctl = np.array([row[colmap[s]] for s in CTRL_SAMPLES], dtype=float)
    recs.append(dict(row=ri, protein_name=name, uniprot=acc_field, uniprot_primary=primary,
                     uniprot_members='|'.join(members), peptide_count=peptide_count,
                     dkd=dkd, ctl=ctl))
wb.close()
print('proteins parsed:', len(recs))

# gene mapping
with open(os.path.join(WORK, 'uniprot_mapping.json'), encoding='utf-8') as f:
    umap = json.load(f)

def gene_of(rec):
    # prefer primary accession gene; fall back to any group member with a gene
    for a in [rec['uniprot_primary']] + rec['uniprot_members'].split('|'):
        if a in umap and umap[a].get('gene'):
            return umap[a]['gene'], a
    return '', rec['uniprot_primary']

stats = []
n_exact = 0
for r in recs:
    dkd, ctl = r['dkd'], r['ctl']
    d_ok, c_ok = dkd[~np.isnan(dkd)], ctl[~np.isnan(ctl)]
    detected_n = int(len(d_ok) + len(c_ok))
    gene, gacc = gene_of(r)
    if len(d_ok) >= 2 and len(c_ok) >= 2:
        med_d, med_c = float(np.median(d_ok)), float(np.median(c_ok))
        l2fc = med_d - med_c
        try:
            u = mannwhitneyu(d_ok, c_ok, alternative='two-sided', method='auto')
            p, used = float(u.pvalue), 'exact'
        except Exception:
            u = mannwhitneyu(d_ok, c_ok, alternative='two-sided')
            p, used = float(u.pvalue), 'asymptotic'
        if used == 'exact':
            n_exact += 1
    else:
        med_d = med_c = l2fc = p = np.nan
        used = ''
    stats.append(dict(gene=gene, uniprot=r['uniprot_primary'], protein_name=r['protein_name'],
                      median_DKD=med_d, median_ctrl=med_c, log2FC=l2fc, wilcox_p=p,
                      test=used, detected_n=detected_n, peptide_count=r['peptide_count']))

df = pd.DataFrame(stats)
print('rows with full stats:', df.log2FC.notna().sum(), '| exact tests:', n_exact)
print('missing values across matrix:', sum(r['dkd'].size - (~np.isnan(r['dkd'])).sum() + r['ctl'].size - (~np.isnan(r['ctl'])).sum() for r in recs))

# percentile of log2FC among all quantified proteins (higher = more up in DKD)
valid = df.log2FC.notna()
df['percentile_log2FC'] = np.nan
df.loc[valid, 'percentile_log2FC'] = df.loc[valid, 'log2FC'].rank(pct=True) * 100

# ---- candidate genes: explicit detection status ----
CANDIDATES = {
    'PDK4': 'Q16654', 'PDK2': 'Q15119',
    'PTGS2': 'P35354', 'NFE2L2': 'Q16236',
    'PDK1': 'Q15118', 'PDK3': 'Q15120',
}
candidate_rows = []
for g, acc in CANDIDATES.items():
    hit = df[df.uniprot == acc]
    if len(hit):
        candidate_rows.append(hit.iloc[0].to_dict())
    else:
        candidate_rows.append(dict(gene=g, uniprot=acc, protein_name='NOT DETECTED in PXD041884',
                                   median_DKD=np.nan, median_ctrl=np.nan, log2FC=np.nan,
                                   wilcox_p=np.nan, test='', detected_n=0, peptide_count=np.nan,
                                   percentile_log2FC=np.nan))
cdf = pd.DataFrame(candidate_rows)
print()
print('=== CANDIDATES ===')
print(cdf[['gene', 'uniprot', 'log2FC', 'wilcox_p', 'detected_n', 'percentile_log2FC']].to_string(index=False))

# ---- context genes (pyruvate / PDH metabolism) ----
CONTEXT = ['PDHA1', 'PDHB', 'DLAT', 'DLD', 'PKM', 'PKLR', 'LDHA', 'LDHB', 'SLC2A1', 'SLC2A2', 'ACLY', 'CS', 'IDH2', 'SDHA']
print()
print('=== CONTEXT (pyruvate/PDH/TCA) ===')
ctx = df[df.gene.isin(CONTEXT)].sort_values('gene')
print(ctx[['gene', 'uniprot', 'median_DKD', 'median_ctrl', 'log2FC', 'wilcox_p', 'detected_n', 'percentile_log2FC']].to_string(index=False))

# dataset-wide summary
print()
print('=== DATASET SUMMARY ===')
print('proteins quantified:', int(valid.sum()))
print('log2FC distribution: median %.3f, IQR [%.3f, %.3f], range [%.3f, %.3f]' % (
    df.loc[valid,'log2FC'].median(), df.loc[valid,'log2FC'].quantile(.25), df.loc[valid,'log2FC'].quantile(.75),
    df.loc[valid,'log2FC'].min(), df.loc[valid,'log2FC'].max()))
print('proteins with wilcox_p < 0.05:', int((df.wilcox_p < 0.05).sum()))
print('proteins down in DKD (log2FC<0):', int((df.loc[valid,'log2FC'] < 0).sum()), '/', int(valid.sum()))

# ---- assemble final deliverable ----
final = pd.concat([cdf, df], ignore_index=True)
final['status'] = np.where(final.protein_name.astype(str).str.startswith('NOT DETECTED'),
                           'not_detected', 'quantified')
cols = ['gene', 'uniprot', 'protein_name', 'status', 'median_DKD', 'median_ctrl', 'log2FC',
        'wilcox_p', 'test', 'detected_n', 'peptide_count', 'percentile_log2FC']
final = final[cols]
out = os.path.join(WORK, 'TableB1_protein_pxd041884.csv')
final.to_csv(out, index=False, float_format='%.6g')
print()
print('WROTE', out, 'rows:', len(final))

# save context table separately for report convenience
ctx_out = os.path.join(WORK, 'TableB1_context_pdh_pyruvate.csv')
pd.concat([cdf, ctx]).to_csv(ctx_out, index=False, float_format='%.6g')
print('WROTE', ctx_out)

# top extremes for narrative
d2 = df[df.log2FC.notna()].copy()
print()
print('Most decreased in DKD:')
print(d2.nsmallest(8, 'log2FC')[['gene','uniprot','log2FC','wilcox_p','percentile_log2FC']].to_string(index=False))
print('Most increased in DKD:')
print(d2.nlargest(8, 'log2FC')[['gene','uniprot','log2FC','wilcox_p','percentile_log2FC']].to_string(index=False))
