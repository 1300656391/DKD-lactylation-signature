"""B1 analysis: lactate-axis protein panel from the existing 1,662-protein matrix.

Extracts DKD-vs-control stats (log2FC, exact Wilcoxon P, percentile) for the
lactate-lactylation axis panel from TableB1_protein_pxd041884.csv (all stats
were pre-computed by analyze_pxd041884.py on the normalized rollup).
Output: TableB1_lactate_axis_panel.csv
"""
import os
import pandas as pd
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
WORK = BASE
df = pd.read_csv(WORK + r'\TableB1_protein_pxd041884.csv')
print('matrix table rows:', len(df), '| quantified:', (df.status == 'quantified').sum())

# (gene, node_class, note) — panel per specification + the module's own candidates
PANEL = [
    ('LDHA', 'lactate_generation', 'lactate dehydrogenase A'),
    ('LDHB', 'lactate_generation', 'lactate dehydrogenase B'),
    ('SLC16A1', 'lactate_transport', 'MCT1'),
    ('SLC16A3', 'lactate_transport', 'MCT4'),
    ('SLC5A12', 'lactate_transport', 'SMCT2 (Na-coupled monocarboxylate)'),
    ('PDK1', 'pyruvate_gate_PDK', ''),
    ('PDK2', 'pyruvate_gate_PDK', 'module candidate'),
    ('PDK3', 'pyruvate_gate_PDK', ''),
    ('PDK4', 'pyruvate_gate_PDK', 'module candidate'),
    ('PDHA1', 'PDH_complex', ''),
    ('PDHB', 'PDH_complex', ''),
    ('DLAT', 'PDH_complex', ''),
    ('DLD', 'PDH_complex', ''),
    ('MPC1', 'mitochondrial_pyruvate_import', ''),
    ('MPC2', 'mitochondrial_pyruvate_import', ''),
    ('PFKFB3', 'glycolysis', '6-phosphofructo-2-kinase/fructose-2,6-bisphosphatase 3'),
    ('PFKL', 'glycolysis', ''),
    ('PFKM', 'glycolysis', ''),
    ('PFKP', 'glycolysis', ''),
    ('PKM', 'glycolysis', ''),
    ('PKLR', 'glycolysis', ''),
    ('ALDOB', 'glycolysis', 'fructose-bisphosphate aldolase B'),
    ('ACLY', 'acetylCoA_source', 'ATP citrate lyase'),
    ('ACSS2', 'acetylCoA_source', 'acetyl-CoA synthetase 2'),
    ('EP300', 'KAT_lactyltransferase_candidate', 'p300 KAT; lysine lactyltransferase candidate'),
    ('HDAC1', 'deacylase', ''),
    ('HDAC2', 'deacylase', ''),
    ('HDAC3', 'deacylase', ''),
    ('SIRT1', 'deacylase', 'NAD-dependent'),
    ('SIRT2', 'deacylase', 'NAD-dependent'),
    ('SIRT3', 'deacylase', 'mitochondrial NAD-dependent'),
    ('SIRT5', 'deacylase', 'mitochondrial; demalonylase/desuccinylase'),
    ('ME1', 'malic_enzyme_NADPH', ''),
    ('ME2', 'malic_enzyme_NADPH', ''),
    ('KRT18', 'tubule_injury_marker', ''),
    ('HIST1H4B', 'histone_H4', 'H4 family quantified as single sequence-identical entry (P62805; gene-level attribution to HIST1H4B/H4C2 impossible)'),
]

# gene name variants in the mapped table (rollup gene symbols came from UniProt)
ALIAS = {'HIST1H4B': ['HIST1H4B', 'H4C2', 'H4C1;']}

rows = []
for gene, cls, note in PANEL:
    variants = ALIAS.get(gene, [gene])
    hit = df[df.gene.isin(variants) & (df.status == 'quantified')]
    if len(hit) == 0:
        hit = df[df.gene.isin(variants)]  # candidate 'not_detected' rows carry status flag
    if len(hit):
        r = hit.iloc[0]
        rows.append(dict(gene=gene, node_class=cls, note=note,
                         status=r.status, uniprot=r.uniprot,
                         protein_name=r.protein_name,
                         median_DKD=r.median_DKD, median_ctrl=r.median_ctrl,
                         log2FC=r.log2FC, wilcox_p=r.wilcox_p,
                         detected_n=r.detected_n, peptide_count=r.peptide_count,
                         percentile_log2FC=r.percentile_log2FC))
    else:
        rows.append(dict(gene=gene, node_class=cls, note=note,
                         status='absent_from_matrix', uniprot='', protein_name='',
                         median_DKD=np.nan, median_ctrl=np.nan, log2FC=np.nan,
                         wilcox_p=np.nan, detected_n=0, peptide_count=np.nan,
                         percentile_log2FC=np.nan))

out = pd.DataFrame(rows)
out.to_csv(WORK + r'\TableB1_lactate_axis_panel.csv', index=False, float_format='%.6g')
print('WROTE TableB1_lactate_axis_panel.csv rows:', len(out))

pd.set_option('display.width', 250)
print()
print(out[['gene', 'node_class', 'status', 'log2FC', 'wilcox_p', 'detected_n',
           'peptide_count', 'percentile_log2FC']].to_string(index=False))

print()
print('=== detection summary ===')
print('quantified:', (out.status == 'quantified').sum(), '/', len(out))
print('significant (p<0.05):')
sig = out[(out.status == 'quantified') & (out.wilcox_p < 0.05)]
print(sig[['gene', 'node_class', 'log2FC', 'wilcox_p', 'percentile_log2FC']].to_string(index=False))
print()
print('not quantified:', out.loc[out.status != 'quantified', ['gene', 'node_class', 'status']].to_string(index=False))
