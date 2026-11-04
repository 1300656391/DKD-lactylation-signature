"""B1 analysis: acetylome landscape — Fritz_Acetyl_Quantitation_FFPE_Area2.xlsx.

289 acetyl-lysine peptides quantified per-sample (DKD n=5 vs control n=7, acetyl-IP
fraction). Never differentially analyzed before (only searched for 4 candidates).
Extracts per-site DKD-vs-control differences on BOTH deposited normalization scales
(Procal-normalized primary; total-protein-adjusted secondary).
Output: TableB1_acetylome_inventory.csv (all sites + stats + pathway flags)
NOTE: acetyl-lysine != lactyl-lysine; this is an inventory of the quantified
acetylation landscape, used only to ask whether PDH-complex/LDH acetylation
is quantified at all.
"""
import os
import openpyxl, os, json, re
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

BASE = os.path.dirname(os.path.abspath(__file__))
DL = os.path.join(BASE, 'downloads')
WORK = BASE

DKD = ['01_D1', '05_D3', '06_D4', '09_D5', '10_D6']
CTRL = ['03_C1', '04_C2', '07_C3', '08_C4', '12_C6', '13_C7', '14_C9']

with open(WORK + r'\uniprot_mapping.json', encoding='utf-8') as f:
    umap = json.load(f)
acc2gene = {a: v.get('gene', '') for a, v in umap.items()}

# fallback gene names for accessions not in the 1,662-protein mapping
FALLBACK = {
    'P62805': 'H4-26/H4C-group (Histone H4)', 'P84243': 'H3-3B/H3.3 (Histone H3)',
    'P68431': 'H3C1/H3.1 (Histone H3)', 'Q16695': 'H3T (Histone H3.1t)',
    'P68871': 'HBB', 'P69905': 'HBA1/HBA2', 'P02042': 'HBD',
    'Q71DI3': 'H3C2 (Histone H3)', 'P0C0S5': 'H2A group', 'P23527': 'HIST1H1B',
    'Q92522': 'H1FX/H1 group', 'P16403': 'H1F0', 'P07305': 'H1F0',
}

def classify(gene, pname):
    g = gene.upper()
    p = pname.upper()
    sets = {
        'LDH': {'LDHA', 'LDHB', 'LDHAL6A'},
        'PDH_complex': {'PDHA1', 'PDHB', 'DLAT', 'DLD', 'PDHX', 'PDP1', 'PDP2',
                        'PDK1', 'PDK2', 'PDK3', 'PDK4'},
        'glycolysis_gluconeogenesis': {'PKM', 'PKLR', 'ALDOA', 'ALDOB', 'ALDOC', 'GAPDH',
                                       'PGK1', 'PGK2', 'ENO1', 'ENO2', 'ENO3', 'TPI1',
                                       'PFKL', 'PFKM', 'PFKP', 'PFKFB3', 'HK1', 'HK2',
                                       'G6PC', 'FBP1', 'GYS1'},
        'TCA_cycle': {'CS', 'ACO1', 'ACO2', 'IDH1', 'IDH2', 'IDH3A', 'IDH3B', 'IDH3G',
                      'OGDH', 'DLST', 'SUCLA2', 'SUCLG1', 'SUCLG2', 'SDHA', 'SDHB',
                      'FH', 'MDH1', 'MDH2', 'ACL', 'ACLY'},
        'ETC_OXPHOS': {'ATP5F1A', 'ATP5F1B', 'ATP5PB', 'ATP5MC1', 'ATP5PB', 'ATP6V1A',
                       'NDUFS1', 'NDUFV1', 'NDUFS3', 'UQCRC2', 'UQCRC1', 'COX5A',
                       'NDUFA10', 'ATP5F1C', 'ATP5F1D', 'ATP5ME'},
        'ketone_body': {'BDH1', 'BDH2', 'HMGCL', 'HMGCS2', 'ACAT1'},
        'amino_acid metabolism': {'GLUD1', 'GLUL', 'GOT1', 'GOT2', 'OAT', 'ALDH2',
                                  'ALDH4A1', 'PCCB', 'MUT', 'PRODH', 'BCAT2'},
        'fatty_acid': {'ACAA2', 'ACADM', 'ACADL', 'ACADS', 'ACAT1', 'HADHA', 'HADHB',
                       'ETFA', 'ETFB', 'ETFDH', 'CPT2'},
    }
    for cls, gs in sets.items():
        # symbol prefix match (handles ATP5A1 legacy vs ATP5F1A new naming)
        for s in gs:
            if g == s or g.startswith(s + ';'):
                return cls
    # legacy-name catches
    if g.startswith('ATP5') or g.startswith('NDUF') or g.startswith('UQCR') or g.startswith('COX'):
        return 'ETC_OXPHOS'
    if g.startswith('HIST') or g.startswith('H1F') or g.startswith('H2A') or g.startswith('H2B') \
       or g.startswith('H3') or g.startswith('H4') or g.startswith('MACROH2'):
        return 'histone'
    if 'HISTONE' in p and ('H3' in p or 'H4' in p or 'H2A' in p or 'H2B' in p or 'H1' in p):
        return 'histone'
    if g.startswith('HBB') or g.startswith('HBA') or g.startswith('HBD') or g.startswith('HBG'):
        return 'hemoglobin'
    if g.startswith('ALB'):
        return 'albumin'
    return 'other'

fn = os.path.join(DL, 'Fritz_Acetyl_Quantitation_FFPE_Area2.xlsx')
wb = openpyxl.load_workbook(fn, read_only=True, data_only=True)
ws = wb['AcK Peptides Norm Procal ADJ TP']
rows = ws.iter_rows(values_only=True)
hdr = list(next(rows))

def find_col(prefix, sample):
    want = f'{prefix} Fritz_Acetyl_{sample}'
    for i, h in enumerate(hdr):
        hs = str(h).strip()
        if hs == want or hs.startswith(want + ':'):
            return i
    raise KeyError(want)

col_procal = {s: find_col('[Area]', s) for s in DKD + CTRL}
col_tpadj = {s: find_col('[Area] UnNorm TP Adj', s) for s in DKD + CTRL}

recs = []
for r in rows:
    if r[0] is None:
        continue
    name = str(r[0])
    f = name.split('&')
    pname, acc = f[0], f[1]
    seq = f[3] if len(f) > 3 else ''
    sites = f[4] if len(f) > 4 else ''
    gene = acc2gene.get(acc.split('|')[0], '') or FALLBACK.get(acc.split('|')[0], '')
    def vals(cols):
        return np.array([float(r[c]) if r[c] is not None else np.nan for c in cols], dtype=float)
    pc_d, pc_c = vals([col_procal[s] for s in DKD]), vals([col_procal[s] for s in CTRL])
    tp_d, tp_c = vals([col_tpadj[s] for s in DKD]), vals([col_tpadj[s] for s in CTRL])
    recs.append(dict(protein_name=pname, uniprot=acc, gene=gene, sequence=seq, sites=sites,
                     procal_dkd=pc_d, procal_ctrl=pc_c, tp_dkd=tp_d, tp_ctrl=tp_c))
wb.close()
print('acetyl sites parsed:', len(recs))

def stats_block(d, c):
    """log2 stats on detected (>0) values; DKD vs ctrl."""
    ld = np.log2(d[d > 0])
    lc = np.log2(c[c > 0])
    n_d, n_c = len(ld), len(lc)
    if n_d >= 3 and n_c >= 3:
        l2fc = float(np.median(ld) - np.median(lc))
        p = float(mannwhitneyu(ld, lc, alternative='two-sided').pvalue)
    else:
        l2fc = np.nan
        p = np.nan
    return n_d, n_c, l2fc, p

rows_out = []
for rec in recs:
    gene = rec['gene']
    cls = classify(gene, rec['protein_name'])
    nd_pc, nc_pc, l2_pc, p_pc = stats_block(rec['procal_dkd'], rec['procal_ctrl'])
    nd_tp, nc_tp, l2_tp, p_tp = stats_block(rec['tp_dkd'], rec['tp_ctrl'])
    rows_out.append(dict(
        protein_name=rec['protein_name'], uniprot=rec['uniprot'], gene=gene,
        pathway=cls, sequence=rec['sequence'], sites=rec['sites'],
        detected_n_DKD_procal=nd_pc, detected_n_ctrl_procal=nc_pc,
        log2FC_procal=l2_pc, wilcox_p_procal=p_pc,
        detected_n_DKD_tpadj=nd_tp, detected_n_ctrl_tpadj=nc_tp,
        log2FC_tpadj=l2_tp, wilcox_p_tpadj=p_tp,
    ))
tab = pd.DataFrame(rows_out)

# BH within testable sites, per scale
for scale in ['procal', 'tpadj']:
    m = tab[f'wilcox_p_{scale}'].notna()
    k = int(m.sum())
    ps = tab.loc[m, f'wilcox_p_{scale}'].values
    o = np.argsort(ps)
    q = np.minimum.accumulate((ps[o] * k / np.arange(1, k + 1))[::-1])[::-1]
    bh = np.empty(k)
    bh[o] = np.clip(q, 0, 1)
    tab.loc[m, f'p_BH_{scale}'] = bh
print('testable sites (>=3 detected/group): procal', int(tab.wilcox_p_procal.notna().sum()),
      '| tpadj', int(tab.wilcox_p_tpadj.notna().sum()))

tab.to_csv(WORK + r'\TableB1_acetylome_inventory.csv', index=False, float_format='%.6g')
print('WROTE TableB1_acetylome_inventory.csv rows:', len(tab))

# ---------- summary ----------
pd.set_option('display.width', 260)
print()
print('=== pathway coverage of the 289 acetyl sites ===')
print(tab.groupby('pathway').size().sort_values(ascending=False).to_string())

print()
print('=== PDH-complex / LDH acetyl sites (the specific ask) ===')
focus = tab[tab.pathway.isin(['PDH_complex', 'LDH'])]
if len(focus) == 0:
    print('NONE quantified')
else:
    print(focus[['gene', 'protein_name', 'sites', 'sequence', 'pathway',
                 'log2FC_procal', 'wilcox_p_procal', 'log2FC_tpadj', 'wilcox_p_tpadj']].to_string(index=False))

print()
print('=== top 15 DKD-increased acetyl sites (Procal scale, testable) ===')
t = tab[tab.wilcox_p_procal.notna()].sort_values('log2FC_procal', ascending=False)
print(t.head(15)[['gene', 'sites', 'pathway', 'log2FC_procal', 'wilcox_p_procal', 'p_BH_procal',
                  'log2FC_tpadj', 'wilcox_p_tpadj']].to_string(index=False))
print()
print('=== top 15 DKD-decreased acetyl sites (Procal scale, testable) ===')
print(t.tail(15).iloc[::-1][['gene', 'sites', 'pathway', 'log2FC_procal', 'wilcox_p_procal', 'p_BH_procal',
                             'log2FC_tpadj', 'wilcox_p_tpadj']].to_string(index=False))

print()
print('=== metabolic-enzyme subset summary (glycolysis/TCA/PDH/LDH/ETC/ketone/FA/AA) ===')
met = tab[tab.pathway.isin(['glycolysis_gluconeogenesis', 'TCA_cycle', 'PDH_complex', 'LDH',
                            'ETC_OXPHOS', 'ketone_body', 'fatty_acid', 'amino_acid metabolism'])]
print('metabolic sites:', len(met), '| testable:', int(met.wilcox_p_procal.notna().sum()))
tm = met[met.wilcox_p_procal.notna()]
print('median log2FC (procal): %+.3f | positive: %d/%d | p<0.05: %d | BH<0.05: %d' % (
    tm.log2FC_procal.median(), (tm.log2FC_procal > 0).sum(), len(tm),
    (tm.wilcox_p_procal < 0.05).sum(), (tm.p_BH_procal < 0.05).sum()))

print()
print('=== histone acetyl sites ===')
h = tab[tab.pathway == 'histone']
print(h[['gene', 'protein_name', 'sites', 'log2FC_procal', 'wilcox_p_procal',
         'log2FC_tpadj', 'wilcox_p_tpadj']].to_string(index=False))

print()
print('=== scale concordance (procal vs tpadj log2FC, testable on both) ===')
both = tab[tab.log2FC_procal.notna() & tab.log2FC_tpadj.notna()]
from scipy.stats import spearmanr
print('n =', len(both), 'spearman rho = %.3f' % spearmanr(both.log2FC_procal, both.log2FC_tpadj).statistic,
      '| sign concordant: %d/%d' % ((np.sign(both.log2FC_procal) == np.sign(both.log2FC_tpadj)).sum(), len(both)))

# unmapped genes?
un = tab[tab.gene == '']
print()
print('accessions with no gene mapping:', len(un), un.uniprot.str.split('|').str[0].unique()[:20])

# ---------- site change beyond protein-abundance change ----------
prot = pd.read_csv(WORK + r'\TableB1_protein_pxd041884.csv')
prot = prot[prot.status == 'quantified'].set_index('gene')['log2FC']
print()
print('=== site-level vs protein-level log2FC (does the acetyl change exceed protein loss?) ===')
for g in ['LDHA', 'LDHB', 'DLD', 'DLAT', 'PDHA1', 'PDHB', 'IDH2', 'GLUD1', 'GATM', 'SLC25A4', 'PKM', 'MDH2', 'ALDOA']:
    sub = tab[tab.gene.str.split(';').str[0] == g]
    if len(sub) == 0 or g not in prot.index:
        continue
    for _, r in sub.iterrows():
        print(f'{g:8s} site {r.sites:20s} protein_l2fc={prot[g]:+.3f} | acetyl_l2fc procal={r.log2FC_procal:+.3f} '
              f'tpadj={r.log2FC_tpadj:+.3f} | EXCESS(procal)={r.log2FC_procal - prot[g]:+.3f} '
              f'EXCESS(tpadj)={r.log2FC_tpadj - prot[g]:+.3f}')
