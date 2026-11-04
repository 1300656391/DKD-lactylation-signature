import os
import requests, json, os, csv, time

BASE = os.path.dirname(os.path.abspath(__file__))
WORK = BASE
DL = os.path.join(WORK, 'downloads')
s = requests.Session()
s.headers.update({'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) DKD-study-B3'})
API = 'https://www.ilincs.org/api/ilincsR/findConcordancesSC'

rows = []
with open(os.path.join(BASE, 'ESM_3.csv'), encoding='utf-8-sig') as f:
    for r in csv.DictReader(f):
        try:
            rows.append((r['gene'].strip(), float(r['logFC']), float(r['adj.P.Val'])))
        except Exception:
            pass
ups_all = sorted([r for r in rows if r[1] > 0], key=lambda r: (r[2], -abs(r[1])))
downs_all = sorted([r for r in rows if r[1] < 0], key=lambda r: (r[2], -abs(r[1])))
gu_all = [g for g, v, p in ups_all]
gd_all = [g for g, v, p in downs_all]

def query_save(lib, up, dn, fname, tries=4):
    sp = json.dumps({'genesUp': up, 'genesDown': dn})
    for a in range(tries):
        try:
            r = s.post(API, data={'mode': 'UpDn', 'libName': lib, 'signatureProfile': sp}, timeout=420)
            j = r.json()
            rem = str((j.get('Remark') or [''])[0])
            tab = j.get('concordanceTable') or []
            n = tab[0].get('nGenes') if tab and isinstance(tab[0], dict) else None
            print('%s: HTTP %s remark=%r nGenes=%s nTable=%d' % (fname, r.status_code, rem[:60], n, len(tab)))
            if tab:
                with open(os.path.join(DL, fname), 'w', encoding='utf-8') as f:
                    json.dump(j, f)
                print('  saved %s (%d KB)' % (fname, os.path.getsize(os.path.join(DL, fname)) // 1024))
                return j
        except Exception as e:
            print('  retry %d: %s' % (a + 1, repr(e)[:90]))
            time.sleep(8)
    return None

# Variant A (up extended to 600 to obtain mappable up genes)
query_save('LIB_5', gu_all[:600], gd_all[:150], 'concordances_LIB_5_final_top600up_top150dn.json')

# Variant B: full-DEG robustness
query_save('LIB_5', gu_all, gd_all, 'concordances_LIB_5_final_fullDEG.json')
