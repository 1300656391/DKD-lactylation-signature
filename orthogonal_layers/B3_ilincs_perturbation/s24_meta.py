import os
import requests, json, os, time

BASE = os.path.dirname(os.path.abspath(__file__))
WORK = BASE
DL = os.path.join(WORK, 'downloads')
s = requests.Session()
s.headers.update({'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) DKD-study-B3'})
API = 'https://www.ilincs.org/api'

def ids_from_ranked(fn, key):
    tab = json.load(open(os.path.join(DL, fn), encoding='utf-8'))
    return [r[key] for r in tab]

need = set()
for fn in ['concordances_LIB_5_final_top600up_top150dn_rankedNeg.json',
           'concordances_LIB_5_final_fullDEG_rankedNeg.json']:
    need.update(ids_from_ranked(fn, 'signatureID'))
need.update(ids_from_ranked('concordances_LIB_2_rankedNeg.json', 'signatureID'))
need = sorted(need)
print('total signature IDs needing metadata:', len(need))

# try batch loopback filter
flt = json.dumps({'where': {'signatureid': {'inq': need}}})
try:
    r = s.get(API + '/SignatureMeta', params={'filter': flt}, timeout=180)
    print('batch GET /SignatureMeta ->', r.status_code, len(r.content), 'bytes')
    if r.status_code == 200:
        metas = r.json()
        print('  got', len(metas), 'metas')
        json.dump(metas, open(os.path.join(DL, 'signaturemeta_batch.json'), 'w'))
except Exception as e:
    print('batch ERR', repr(e)[:120])

# fallback: individual
got = {}
if os.path.exists(os.path.join(DL, 'signaturemeta_batch.json')):
    for m in json.load(open(os.path.join(DL, 'signaturemeta_batch.json'), encoding='utf-8')):
        got[m.get('signatureid')] = m
todo = [i for i in need if i not in got]
print('still needing individual fetch:', len(todo))
for i, sid in enumerate(todo):
    for a in range(3):
        try:
            r = s.get(API + '/SignatureMeta/' + sid, timeout=60)
            if r.status_code == 200:
                got[sid] = r.json()
                break
            elif r.status_code == 404:
                got[sid] = {'signatureid': sid, 'missing': True}
                break
        except Exception:
            time.sleep(3)
    if (i + 1) % 20 == 0:
        print('  %d/%d fetched' % (i + 1, len(todo)))
    time.sleep(0.4)

json.dump(got, open(os.path.join(DL, 'signaturemeta_all.json'), 'w'))
ok = sum(1 for v in got.values() if not v.get('missing'))
print('metadata collected: %d ok, %d missing, total %d' % (ok, len(got) - ok, len(need)))
# show sample
for sid in need[:3]:
    m = got.get(sid, {})
    print(sid, '->', {k: m.get(k) for k in ['compound', 'cellline', 'concentration', 'lincspertid', 'libraryid']})
