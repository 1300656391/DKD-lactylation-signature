"""s34_build_table.py - Build deliverable TableB3_ilincs_reversers.csv (no API calls).

Primary (definitive) query = concordances_LIB_5_final_fullDEG.json:
  POST /api/ilincsR/findConcordancesSC  mode=UpDn  libName=LIB_5
  signatureProfile = {genesUp: all 1259 up DEGs, genesDown: all 1386 down DEGs}
  -> 114 genes mapped, 14,337 signatures scored, 9,565 negative.

Table = top 50 reversing perturbations (most-negative similarity) joined to
signaturemeta_all.json; per-compound median similarity computed as in
s28_build_table.py (median across the compound's signatures among the 120
strongest reversing signatures of the definitive query).
"""
import os
import json, os, csv, statistics, re

BASE = os.path.dirname(os.path.abspath(__file__))
WORK = BASE
DL = os.path.join(WORK, 'downloads')

metas = json.load(open(os.path.join(DL, 'signaturemeta_all.json'), encoding='utf-8'))
rnB = json.load(open(os.path.join(DL, 'concordances_LIB_5_final_fullDEG_rankedNeg.json'), encoding='utf-8'))
rnA = json.load(open(os.path.join(DL, 'concordances_LIB_5_final_top600up_top150dn_rankedNeg.json'), encoding='utf-8'))
idsA = {r['signatureID'] for r in rnA}

def parse_lincs_sid(lsid):
    m = re.match(r'^([^:]+)_(.+?)_(\d+[HM]):([^:]+):(.+)$', lsid or '')
    if m:
        return m.group(2), m.group(5), m.group(3)   # cellline, dose, time
    return '', '', ''

# per-compound aggregation over the definitive top-120 set (s28 semantics)
by_compound = {}
for r in rnB:
    m = metas.get(r['signatureID'], {})
    p = parse_lincs_sid(r.get('lincsSigID'))
    comp = m.get('compound') or p[1] or r.get('lincsSigID') or r['signatureID']
    by_compound.setdefault(comp, []).append(r['similarity'])
comp_stats = {c: {'n': len(s), 'median': statistics.median(s)} for c, s in by_compound.items()}

out = os.path.join(WORK, 'TableB3_ilincs_reversers.csv')
with open(out, 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f)
    w.writerow(['rank', 'compound_name', 'lincsSigID', 'cellline', 'dose', 'timepoint',
                'similarity', 'pValue', 'nGenes_mapped', 'MoA_category',
                'n_signatures_of_compound_in_top120', 'compound_median_similarity_top120',
                'also_in_top120_variantA'])
    for i, r in enumerate(rnB[:50], 1):
        m = metas.get(r['signatureID'], {})
        cell, dose, tp = parse_lincs_sid(r.get('lincsSigID'))
        comp = m.get('compound') or dose or r.get('lincsSigID')
        st = comp_stats.get(comp, {'n': 1, 'median': r['similarity']})
        w.writerow([i, comp, r.get('lincsSigID') or '',
                    m.get('cellline') or cell,
                    m.get('concentration') or dose, m.get('time') or tp,
                    '%.4f' % r['similarity'], '%.3g' % (r.get('pValue') or float('nan')),
                    r.get('nGenes'), m.get('integratedMoas') or '',
                    st['n'], '%.4f' % st['median'],
                    'yes' if r['signatureID'] in idsA else 'no'])
print('wrote', out)
for i, r in enumerate(rnB[:50], 1):
    m = metas.get(r['signatureID'], {})
    print('%2d sim=%+.4f p=%-8.2g %-40s %-8s %-10s %-5s %s' % (
        i, r['similarity'], r.get('pValue') or -1,
        str(m.get('compound'))[:40], str(m.get('cellline'))[:8],
        str(m.get('concentration'))[:10], str(m.get('time'))[:5],
        (str(m.get('integratedMoas')) or '')[:44]))
