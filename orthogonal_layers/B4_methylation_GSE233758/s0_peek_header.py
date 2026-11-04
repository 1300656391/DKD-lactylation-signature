import os
import gzip

BASE = os.path.dirname(os.path.abspath(__file__))
p = os.path.join(BASE, 'downloads', 'GSE233758_series_matrix.txt.gz')
out = os.path.join(BASE, 'header_meta.txt')
n = 0
with gzip.open(p, 'rt', encoding='utf-8', errors='replace') as f, open(out, 'w', encoding='utf-8') as w:
    for line in f:
        if line.startswith('!series_matrix_table_begin'):
            break
        w.write(line)
        n += 1
print('header lines written:', n)
