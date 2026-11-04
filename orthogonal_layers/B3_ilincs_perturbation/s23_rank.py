import os
import json, os, statistics

BASE = os.path.dirname(os.path.abspath(__file__))
DL = os.path.join(BASE, 'downloads')

for fn in ['concordances_LIB_5_final_top600up_top150dn.json',
           'concordances_LIB_5_final_fullDEG.json',
           'concordances_LIB_2.json']:
    j = json.load(open(os.path.join(DL, fn), encoding='utf-8'))
    tab = j['concordanceTable']
    tab = [r for r in tab if isinstance(r, dict) and 'similarity' in r]
    sims = [r['similarity'] for r in tab]
    neg = sorted([r for r in tab if r['similarity'] < 0], key=lambda r: r['similarity'])
    print('== %s: n=%d, neg=%d, sim range %.3f..%.3f' % (fn, len(tab), len(neg), min(sims), max(sims)))
    for r in neg[:10]:
        print('   sim=%+.4f p=%-10.3g nG=%-4s %s' % (r['similarity'], r.get('pValue', -1), str(r.get('nGenes')), r.get('lincsSigID') or r.get('signatureID')))
    # save ranked negative list top 120 for metadata fetch
    out = fn.replace('.json', '_rankedNeg.json')
    json.dump(neg[:120], open(os.path.join(DL, out), 'w'))
    print('   saved', out)
