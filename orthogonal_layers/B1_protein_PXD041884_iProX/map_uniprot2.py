"""Map accessions to gene symbols via direct UniProt search API (no job polling)."""
import os
import openpyxl, os, json, time, urllib.request, urllib.parse

BASE = os.path.dirname(os.path.abspath(__file__))
DL = os.path.join(BASE, 'downloads')
WORK = BASE

wb = openpyxl.load_workbook(os.path.join(DL, 'Secocnd_Extraction_GPQ_Norm_TS_All_Peptides_FFPE_Protein_Rollup_Workbook.xlsx'), read_only=True, data_only=True)
ws = wb['All Protien Rolled Up']
rows = ws.iter_rows(values_only=True)
next(rows)
all_accs = set()
for row in rows:
    acc = str(row[1]).strip() if row[1] else ''
    for a in acc.split('|'):
        all_accs.add(a)
wb.close()
acc_list = sorted(all_accs)
print('unique accessions:', len(acc_list), flush=True)

results = {}
CHUNK = 40
for i in range(0, len(acc_list), CHUNK):
    chunk = acc_list[i:i+CHUNK]
    q = ' OR '.join('accession:' + a for a in chunk)
    url = 'https://rest.uniprot.org/uniprotkb/search?' + urllib.parse.urlencode(
        {'query': q, 'format': 'tsv', 'fields': 'accession,gene_primary,protein_name', 'size': '500'})
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={'Accept': 'text/tsv'})
            with urllib.request.urlopen(req, timeout=90) as r:
                tsv = r.read().decode()
            lines = tsv.splitlines()
            cols = lines[0].split('\t')
            for l in lines[1:]:
                parts = (l + '\t').split('\t')
                d = dict(zip(cols, parts))
                acc = d.get('Entry', '')
                genes = d.get('Gene Names (primary)', '') or ''
                prim = genes.split()[0] if genes else ''
                results[acc] = {'gene': prim, 'all_genes': genes, 'name': d.get('Protein names', '')}
            print(f'chunk {i//CHUNK}: +{len(lines)-1} entries (cum {len(results)})', flush=True)
            break
        except Exception as e:
            print('chunk', i//CHUNK, 'attempt', attempt, 'ERR', repr(e)[:150], flush=True)
            time.sleep(5)
    time.sleep(1)

with open(os.path.join(WORK, 'uniprot_mapping.json'), 'w', encoding='utf-8') as f:
    json.dump(results, f, ensure_ascii=False, indent=1)
print('total mapped accessions:', len(results), flush=True)
unmapped = [a for a in acc_list if a not in results]
print('unmapped:', len(unmapped), unmapped[:20], flush=True)
for k in ['P14618', 'P08559', 'P11177', 'P10515', 'P04446', 'Q16654', 'Q15119']:
    print(k, '->', results.get(k, {}).get('gene'), '|', str(results.get(k, {}).get('name'))[:60], flush=True)
