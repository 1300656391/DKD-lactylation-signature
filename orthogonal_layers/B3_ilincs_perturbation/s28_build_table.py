import os
import json, os, csv, re, statistics

BASE = os.path.dirname(os.path.abspath(__file__))
WORK = BASE
DL = os.path.join(WORK, 'downloads')

metas = json.load(open(os.path.join(DL, 'signaturemeta_all.json'), encoding='utf-8'))

def load_ranked(fn):
    return json.load(open(os.path.join(DL, fn), encoding='utf-8'))

def parse_lincs_sid(lsid):
    """e.g. CPC006_A673_6H:BRD-A36275421-001-11-1:10 -> dataset, cell, time, perturbagen, dose"""
    m = re.match(r'^([^:]+)_(.+?)_(\d+[HM]):([^:]+):(.+)$', lsid or '')
    if m:
        return {'dataset': m.group(1), 'cellline': m.group(2), 'time': m.group(3),
                'perturbagen': m.group(4), 'dose': m.group(5)}
    return {'dataset': '', 'cellline': '', 'time': '', 'perturbagen': '', 'dose': ''}

def enrich(rows, lib):
    out = []
    for r in rows:
        sid = r.get('signatureID')
        m = metas.get(sid, {})
        p = parse_lincs_sid(r.get('lincsSigID'))
        out.append({
            'library': lib,
            'signatureID': sid,
            'lincsSigID': r.get('lincsSigID') or '',
            'similarity': r['similarity'],
            'pValue': r.get('pValue'),
            'nGenes': r.get('nGenes'),
            'compound': m.get('compound') or p['perturbagen'],
            'cellline': m.get('cellline') or p['cellline'],
            'tissue': m.get('tissue'),
            'concentration': m.get('concentration') or p['dose'],
            'timepoint': (m.get('time') or p['time'] or ''),
            'moa': m.get('integratedMoas') or '',
            'lincspertid': m.get('lincspertid') or '',
            'pubChemID': m.get('pubChemID') or '',
        })
    return out

lib5_A = enrich(load_ranked('concordances_LIB_5_final_top600up_top150dn_rankedNeg.json'), 'LIB_5')
lib5_B = enrich(load_ranked('concordances_LIB_5_final_fullDEG_rankedNeg.json'), 'LIB_5')
lib2 = enrich(load_ranked('concordances_LIB_2_rankedNeg.json'), 'LIB_2')

print('LIB_5 variant A (top600up/top150dn) neg rows:', len(lib5_A))
print('LIB_5 variant B (full DEG)          neg rows:', len(lib5_B))
print('LIB_2 neg rows:', len(lib2))

# ---- per-compound aggregation for LIB_5 variant A (primary) ----
def agg(rows):
    by_compound = {}
    for r in rows:
        key = r['compound']
        by_compound.setdefault(key, []).append(r)
    agg_rows = []
    for comp, rs in by_compound.items():
        sims = [x['similarity'] for x in rs]
        best = min(rs, key=lambda x: x['similarity'])
        agg_rows.append({
            'compound': comp,
            'n_signatures': len(rs),
            'best_similarity': min(sims),
            'best_pValue': best['pValue'],
            'best_signatureID': best['signatureID'],
            'best_lincsSigID': best['lincsSigID'],
            'best_cellline': best['cellline'],
            'best_concentration': best['concentration'],
            'best_timepoint': best['timepoint'],
            'median_similarity': statistics.median(sims),
            'moa': best['moa'],
            'moa_any': '|'.join(sorted({x['moa'] for x in rs if x['moa']})),
            'lincspertid': best['lincspertid'],
        })
    agg_rows.sort(key=lambda r: r['best_similarity'])
    return agg_rows

aggA = agg(lib5_A)
aggB = agg(lib5_B)
agg2 = agg(lib2)

print('\n== LIB_5 variant A: top 25 compounds by best (most negative) similarity ==')
for a in aggA[:25]:
    print('  %-42s n=%-3d best=%+.4f med=%+.4f p=%-9.3g moa=%s' % (a['compound'][:42], a['n_signatures'], a['best_similarity'], a['median_similarity'], a['best_pValue'] or -1, (a['moa'] or '')[:60]))

print('\n== LIB_5 variant B (full DEG): top 15 ==')
for a in aggB[:15]:
    print('  %-42s n=%-3d best=%+.4f med=%+.4f moa=%s' % (a['compound'][:42], a['n_signatures'], a['best_similarity'], a['median_similarity'], (a['moa'] or '')[:60]))

print('\n== LIB_2: top 15 ==')
for a in agg2[:15]:
    print('  %-42s n=%-3d best=%+.4f med=%+.4f moa=%s' % (a['compound'][:42], a['n_signatures'], a['best_similarity'], a['median_similarity'], (a['moa'] or '')[:60]))

# variant agreement on compounds
topA = {a['compound'] for a in aggA[:50]}
topB = {a['compound'] for a in aggB[:50]}
print('\nvariant A/B top-50 compound overlap:', len(topA & topB))

# save enriched rows + aggregates
json.dump({'lib5A': lib5_A, 'lib5B': lib5_B, 'lib2': lib2}, open(os.path.join(WORK, 'enriched_rows.json'), 'w'))
json.dump({'aggA': aggA, 'aggB': aggB, 'agg2': agg2}, open(os.path.join(WORK, 'aggregate_compounds.json'), 'w'), indent=1)
print('saved enriched_rows.json and aggregate_compounds.json')
