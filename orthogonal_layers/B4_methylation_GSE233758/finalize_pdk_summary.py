import os
# -*- coding: utf-8 -*-
# Restructure the SUMMARY rows of TableB4_pdk_methylation.csv into the self-documenting
# best-probe scheme (run after s5_make_table.py):
#   - probe rows (cgID starts with 'cg'): unchanged; best_probe_* columns empty.
#   - SUMMARY rows: p = region-level t-test p (all probes, unadjusted; region_level_stats.csv);
#     p_BH = BH across the four region-level tests;
#     best_probe_id/p/p_BH = single best (min-p) probe of the category;
#     best_probe_p_batch_adj = minimum batch-adjusted p among the category's probes
#     (for two PDK4 categories this comes from a different probe than best_probe_id).
# All existing rows and values are preserved; nothing deleted.

BASE = os.path.dirname(os.path.abspath(__file__))
import sys, io, csv
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
import numpy as np

B4 = os.path.join(BASE, 'TableB4_pdk_methylation.csv')
RLS = os.path.join(BASE, 'region_level_stats.csv')

rows = list(csv.DictReader(open(B4, encoding='utf-8')))
fieldnames = list(rows[0].keys())
region = list(csv.DictReader(open(RLS, encoding='utf-8')))
rmap = {(r['gene'], r['region'], r['set']): r for r in region}

# BH across the four region-level 'all' tests
keys = [k for k in rmap if k[2] == 'all']
p_all = sorted((float(rmap[k]['p']), k) for k in keys)
m = len(p_all)
prev = 1.0
q_region = {}
for i in range(m - 1, -1, -1):
    prev = min(prev, p_all[i][0] * m / (i + 1))
    q_region[p_all[i][1]] = prev
print('region-level BH across 4 tests:')
for k, q in sorted(q_region.items(), key=lambda kv: kv[1]):
    print(f'  {k}: p={float(rmap[k]["p"]):.6g} -> q={q:.6g}')

probes = [r for r in rows if r['cgID'].startswith('cg')]
summ = [r for r in rows if not r['cgID'].startswith('cg')]

# capture the previous SUMMARY values and the min-p / min-batch-adj probe ids
log = []
for s in summ:
    sub = [r for r in probes if r['gene'] == s['gene'] and r['region_category'] == s['region_category']]
    best = min(sub, key=lambda r: float(r['p']))
    bestb = min(sub, key=lambda r: float(r['p_batch_adj']))
    log.append((s['cgID'], s['p'], s['p_BH'], s['p_batch_adj'], best['cgID'], bestb['cgID']))

newrows = []
for r in rows:
    nr = dict(r)
    if r['cgID'].startswith('cg'):
        nr['best_probe_id'] = ''
        nr['best_probe_p'] = ''
        nr['best_probe_p_BH'] = ''
        nr['best_probe_p_batch_adj'] = ''
    else:
        k = (r['gene'], r['region_category'], 'all')
        sub = [q for q in probes if q['gene'] == r['gene'] and q['region_category'] == r['region_category']]
        best = min(sub, key=lambda q: float(q['p']))
        nr['best_probe_id'] = best['cgID']
        nr['best_probe_p'] = r['p']            # old value preserved
        nr['best_probe_p_BH'] = r['p_BH']      # old value preserved
        nr['best_probe_p_batch_adj'] = r['p_batch_adj']  # old value preserved
        nr['p'] = rmap[k]['p']
        nr['p_BH'] = repr(q_region[k])
        nr['p_batch_adj'] = rmap[k]['p_batch_adj']
    newrows.append(nr)

newfields = []
for f in fieldnames:
    newfields.append(f)
    if f == 'p_BH':
        newfields += ['best_probe_id', 'best_probe_p', 'best_probe_p_BH']
    if f == 'p_batch_adj':
        newfields.append('best_probe_p_batch_adj')

with open(B4, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=newfields)
    w.writeheader()
    w.writerows(newrows)

print('\nprevious SUMMARY values preserved into best_probe_* columns:')
for cgID, p, pb, pba, bpid, bpbid in log:
    print(f'  {cgID}: old p={p}, p_BH={pb}, p_batch_adj={pba}; best-probe={bpid}; min-batch-adj-probe={bpbid}')

# verify round-trip
chk = list(csv.DictReader(open(B4, encoding='utf-8')))
print(f'\nre-read: {len(chk)} rows (expect 65), {len(chk[0])} columns (expect 20)')
print('new header:', list(chk[0].keys()))
for r in chk:
    if not r['cgID'].startswith('cg'):
        print(' ', r['cgID'], '| p=', r['p'], '| p_BH=', r['p_BH'], '| p_batch_adj=', r['p_batch_adj'],
              '| best:', r['best_probe_id'], r['best_probe_p'], r['best_probe_p_BH'], r['best_probe_p_batch_adj'])
