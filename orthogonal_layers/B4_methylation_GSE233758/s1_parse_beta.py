"""Parse GSE233758 series matrix -> beta matrix pickle + sample sheet."""
import os
import gzip
import re
import time
import pandas as pd
import numpy as np

BASE = BASE
SM = BASE + r'\downloads\GSE233758_series_matrix.txt.gz'

t0 = time.time()
samples = []
titles = []
in_table = False
header_done = False
rows = []
colnames = []

def split_tabs(line):
    return line.rstrip('\n').split('\t')

with gzip.open(SM, 'rt', encoding='utf-8', errors='replace') as f:
    for line in f:
        if not header_done:
            if line.startswith('!Sample_geo_accession'):
                samples = [x.strip('"') for x in split_tabs(line)[1:]]
            elif line.startswith('!Sample_title'):
                titles = [x.strip('"') for x in split_tabs(line)[1:]]
            elif line.startswith('!Sample_submission_date'):
                subdates = [x.strip('"') for x in split_tabs(line)[1:]]
            elif line.startswith('!series_matrix_table_begin'):
                in_table = True
                header_done = True
                continue
        elif in_table:
            if line.startswith('!series_matrix_table_end'):
                break
            parts = split_tabs(line)
            if parts[0] == '"ID_REF"':
                colnames = [x.strip('"') for x in parts[1:]]
                continue
            rows.append(parts)

print('parse done', time.time() - t0, 's; rows read:', len(rows))
assert len(colnames) == len(samples) == len(titles), (len(colnames), len(samples), len(titles))

# probe ids and value matrix
ids = [r[0].strip('"') for r in rows]
mat = np.empty((len(rows), len(colnames)), dtype=np.float64)
na_count = 0
for i, r in enumerate(rows):
    for j in range(len(colnames)):
        v = r[j + 1].strip('"')
        if v == 'NA' or v == '' or v == 'null':
            mat[i, j] = np.nan
            na_count += 1
        else:
            mat[i, j] = float(v)
print('NA cells:', na_count)

beta = pd.DataFrame(mat, index=ids, columns=samples)
beta.to_pickle(BASE + r'\beta_matrix.pkl')

# sample sheet
source = []
for t in titles:
    source.append('DN' if 'diabetic nephropathy' in t else 'Control')
batch = ['B1' if s.startswith('GSM7434') else 'B2' for s in samples]
ss = pd.DataFrame({
    'gsm': samples,
    'title': titles,
    'group': source,
    'batch': batch,
    'submission_date': subdates,
})
ss.to_csv(BASE + r'\sample_sheet.csv', index=False)
print(ss.groupby(['group', 'batch']).size())
print('shape:', beta.shape)
print('beta global range:', np.nanmin(beta.values), np.nanmax(beta.values))
# per-sample NA and mean beta QC
qc = pd.DataFrame({'n_na': beta.isna().sum(), 'mean_beta': beta.mean(axis=0, skipna=True)})
qc.to_csv(BASE + r'\qc_per_sample.csv')
print(qc.to_string())
