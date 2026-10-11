# -*- coding: utf-8 -*-
# R11-F1: repair results/immune_gene_corr.csv gene column (row indices -> gene symbols).
# Order = A_candidates.csv gene order = the order 09_immune.R looped (intersect preserves it).
# Numeric matrix untouched.
import csv

rows = list(csv.reader(open(r"results\immune_gene_corr.csv", newline="")))
hdr = rows[0]
gi = hdr.index("gene")
cand = [r[0] for r in csv.reader(open(r"results\A_candidates.csv", newline=""))][1:]
assert len(cand) == 23 and len(rows) - 1 == 23, (len(cand), len(rows) - 1)

old = [r[gi] for r in rows[1:]]
for i, r in enumerate(rows[1:]):
    r[gi] = cand[i]

with open(r"results\immune_gene_corr.csv", "w", newline="") as f:
    csv.writer(f).writerows(rows)

print("replaced col idx", gi)
print("old:", old[:5], "...", old[-3:])
print("new:", [r[gi] for r in rows[1:6]], "...", [r[gi] for r in rows[-4:-1]])
print("n cells cols:", len(hdr) - 1)
