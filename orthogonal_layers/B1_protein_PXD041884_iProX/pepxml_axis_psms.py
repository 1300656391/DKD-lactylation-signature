"""B1 supplementary analysis: sub-threshold PSM sweep for NON-quantified lactate-axis genes in the 7 pep.xml files.

Extends the PDK4-style sub-threshold check to every axis-panel gene absent from the
1,662-protein quantitative matrix: did the deposited search-engine output contain ANY
peptide-spectrum match (any score, incl. below the SM>=8 library cut-off)?
"""
import os, re, glob, json, time, urllib.request, urllib.parse
from collections import defaultdict

BASE = os.path.dirname(os.path.abspath(__file__))
DL = os.path.join(BASE, 'downloads')
WORK = BASE

# non-quantified axis genes (from TableB1_lactate_axis_panel status != quantified, minus PDK family already known)
GENES = ['SLC16A1', 'SLC16A3', 'SLC5A12', 'MPC1', 'MPC2', 'PFKFB3', 'PFKP', 'ACLY', 'ACSS2',
         'EP300', 'HDAC1', 'HDAC2', 'HDAC3', 'SIRT1', 'SIRT2', 'SIRT3', 'SIRT5',
         'ME1', 'ME2', 'HIST1H4B', 'PDK1', 'PDK2', 'PDK3', 'PDK4']

# resolve canonical human reviewed accessions via UniProt REST
acc2gene = {}
def query(gene):
    url = ('https://rest.uniprot.org/uniprotkb/search?query=' +
           urllib.parse.quote(f'gene_exact:{gene} AND organism_id:9606 AND reviewed:true') +
           '&fields=accession,gene_primary&format=tsv&size=5')
    try:
        with urllib.request.urlopen(url, timeout=45) as r:
            txt = r.read().decode()
    except Exception as e:
        print('  query failed', gene, e)
        return []
    out = []
    for line in txt.splitlines()[1:]:
        f = line.split('\t')
        if len(f) >= 2:
            out.append((f[0], f[1]))
    return out

print('--- UniProt resolution ---')
for g in GENES:
    hits = query(g)
    time.sleep(0.3)
    if not hits:  # legacy symbols (HIST1H4B): try direct-known accession
        hits = [('P62807', 'H4C2')] if g == 'HIST1H4B' else []
    if hits:
        acc2gene[hits[0][0]] = g
        print(f'{g}: {hits[0][0]} ({hits[0][1]})' + (f'  [+{len(hits)-1} more]' if len(hits) > 1 else ''))
    else:
        print(f'{g}: NO REVIEWED HUMAN ENTRY FOUND')
print('resolved:', len(acc2gene), 'genes')

# parse pep.xml
ATTR = re.compile(r'(\w+)="([^"]*)"')
files = sorted(glob.glob(os.path.join(DL, '*.pep.xml')))
psm = defaultdict(int)               # (gene, file) -> count
best = {}                            # gene -> (score, file, spectrum, peptide)
for fp in files:
    base = os.path.basename(fp)
    with open(fp, encoding='utf-8', errors='replace') as f:
        text = f.read()
    n_sq = 0
    for sq_m in re.finditer(r'<spectrum_query ([^>]*)>(.*?)</spectrum_query>', text, re.S):
        body = sq_m.group(2)
        n_sq += 1
        for hit_m in re.finditer(r'<search_hit ([^>]*?)(/>|>)', body):
            ha = dict(ATTR.findall(hit_m.group(1)))
            prots = [p for p in ha.get('protein', '').split(',') if p and not p.startswith('rev_')]
            genes = {acc2gene[p] for p in prots if p in acc2gene}
            if not genes:
                continue
            hit_full = body[hit_m.start():hit_m.start() + 3000]
            scores = dict(re.findall(r'<search_score name="([^"]*)" value="([^"]*)"', hit_full))
            try:
                sm = float(scores.get('SMscore', scores.get('xcorr', 'nan')))
            except ValueError:
                sm = float('nan')
            for g in genes:
                psm[(g, base)] += 1
                old = best.get(g)
                better = old is None or (sm == sm and (old[0] != old[0] or sm > old[0]))
                if better:
                    best[g] = (sm, base, ha.get('peptide', ''),
                               ha.get('num_matched_ions', ''), scores)
    print(f'{base}: {n_sq} spectrum queries parsed')

print()
print('=== axis-gene PSM evidence in deposited search output (genes NOT in the quantitative matrix) ===')
rows = []
for g in GENES:
    tot = sum(v for (gg, f), v in psm.items() if gg == g)
    if tot:
        sm, f, pep, ions, scores = best[g]
        print(f'{g}: {tot} PSM(s) | best SMscore {sm} ({f}) peptide {pep} ions {ions} '
              f'| scores {dict(list(scores.items())[:4])}')
        rows.append((g, tot, sm, f, pep))
    else:
        print(f'{g}: zero PSMs at any score (searched-but-not-detected)')
json.dump({g: dict(n_psms=sum(v for (gg, f), v in psm.items() if gg == g),
                   best_sm=best.get(g, (None,))[0],
                   best_file=best.get(g, (None, None))[1] if g in best else None,
                   best_peptide=best.get(g, (None, None, ''))[2] if g in best else None)
           for g in GENES},
          open(WORK + r'\pepxml_axis_psms.json', 'w'), indent=1)
print('WROTE pepxml_axis_psms.json')
