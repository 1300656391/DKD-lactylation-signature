import os
import csv, os, json

BASE = os.path.dirname(os.path.abspath(__file__))
SUP = os.path.join(BASE, 'ESM_3.csv')  # full DEG list (2645 rows)
WORK = BASE

rows = []
with open(SUP, encoding='utf-8-sig') as f:
    rd = csv.DictReader(f)
    print('columns:', rd.fieldnames)
    for r in rd:
        try:
            rows.append({
                'gene': r['gene'].strip(),
                'logFC': float(r['logFC']),
                'adj.P.Val': float(r['adj.P.Val']),
            })
        except (ValueError, KeyError, TypeError):
            print('skip row:', r)

print('total DEG rows parsed:', len(rows))

ups = [r for r in rows if r['logFC'] > 0]
downs = [r for r in rows if r['logFC'] < 0]
print('n up (logFC>0):', len(ups), ' n down (logFC<0):', len(downs), ' n zero:', len(rows) - len(ups) - len(downs))

def sort_key(r):
    return (r['adj.P.Val'], -abs(r['logFC']))

ups.sort(key=sort_key)
downs.sort(key=sort_key)

top_up = ups[:150]
top_down = downs[:150]

for tag, data in [('up', top_up), ('down', top_down)]:
    out = os.path.join(WORK, 'B3_signature_%s.csv' % tag)
    with open(out, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['rank', 'gene', 'logFC', 'adj.P.Val', 'direction'])
        for i, r in enumerate(data, 1):
            w.writerow([i, r['gene'], '%.6g' % r['logFC'], '%.6g' % r['adj.P.Val'], ('UP' if tag == 'up' else 'DOWN')])
    print('wrote', out, len(data), 'genes; adjP range %.3g - %.3g' % (data[0]['adj.P.Val'], data[-1]['adj.P.Val']))

# overlap with known DN markers
markers = ['PDK4', 'PDK2', 'SPP1', 'VCAN', 'COL1A1', 'COL3A1', 'FN1', 'LCP1', 'CD86', 'C1QA', 'C1QB', 'FOS', 'DUSP1', 'JUN', 'NR4A1', 'KNG1', 'CCL2', 'TIMP1', 'SERPINE1', 'ITGB2', 'TYROBP', 'FCER1G', 'CTSS']
for m in markers:
    rank_up = next((i for i, r in enumerate(ups, 1) if r['gene'] == m), None)
    rank_dn = next((i for i, r in enumerate(downs, 1) if r['gene'] == m), None)
    hit = next((r for r in rows if r['gene'] == m), None)
    if hit:
        loc = 'UP-list' if rank_up else ('DOWN-list' if rank_dn else 'not in lists')
        pos = rank_up if rank_up else (rank_dn if rank_dn else None)
        in_top = 'IN top-150 %s' % ('UP' if rank_up else 'DOWN') if (rank_up and rank_up <= 150) or (rank_dn and rank_dn <= 150) else 'outside top150'
        print('%-8s logFC=%+.3f adjP=%.3g  %s (rank %s) -> %s' % (m, hit['logFC'], hit['adj.P.Val'], loc, pos, in_top))
    else:
        print('%-8s NOT in DEG table' % m)

# also dump full ranked lists of PDK4 specifically
for m in ['PDK4']:
    for i, r in enumerate(downs, 1):
        if r['gene'] == m:
            print('PDK4 DOWN significance-rank: %d / %d; in top-150: %s' % (i, len(downs), i <= 150))
    for i, r in enumerate(ups, 1):
        if r['gene'] == m:
            print('PDK4 UP significance-rank: %d / %d; in top-150: %s' % (i, len(ups), i <= 150))

# save JSON payload for API query
payload = {
    'genesUp': [r['gene'] for r in top_up],
    'genesDown': [r['gene'] for r in top_down],
}
with open(os.path.join(WORK, 'signature_payload.json'), 'w') as f:
    json.dump(payload, f)
print('payload saved: %d up / %d down' % (len(payload['genesUp']), len(payload['genesDown'])))
