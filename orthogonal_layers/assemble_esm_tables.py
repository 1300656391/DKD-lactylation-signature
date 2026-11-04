# -*- coding: utf-8 -*-
# Assemble the submission supplementary tables ESM_28-ESM_33 from the module output
# CSVs (run the three module pipelines first; see README). Cell text is preserved
# verbatim; tables are merged by column-name union in first-seen order.
import csv
import os

BASE = os.path.dirname(os.path.abspath(__file__))
B1 = os.path.join(BASE, "B1_protein_PXD041884_iProX")
B4 = os.path.join(BASE, "B4_methylation_GSE233758")
B3 = os.path.join(BASE, "B3_ilincs_perturbation")


def read_csv(folder, name):
    with open(os.path.join(folder, name), encoding="utf-8-sig", newline="") as f:
        rows = list(csv.reader(f))
    return rows[0], rows[1:]


def write_csv(path, header, rows):
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f, lineterminator="\r\n")
        w.writerow(header)
        w.writerows(rows)


def concat_by_name(frames):
    """frames: list of (header, rows, extra_dict). Column union in first-seen order;
    extra (e.g. block/cohort) columns are appended after all source columns."""
    header = []
    for h, _, _ in frames:
        for c in h:
            if c not in header:
                header.append(c)
    for _, _, extra in frames:
        for c in extra:
            if c not in header:
                header.append(c)
    out = []
    for h, rows, extra in frames:
        idx = {c: i for i, c in enumerate(h)}
        for r in rows:
            row = [r[idx[c]] if c in idx else "" for c in header]
            for k, v in extra.items():
                row[header.index(k)] = v
            out.append(row)
    return header, out


def as_dict(header, rows):
    return [{h: r[i] for i, h in enumerate(header)} for r in rows]


# ---- ESM_28: FFPE proteome rollup + lactate-axis panel flags + PDH/pyruvate context ----
h_roll, r_roll = read_csv(B1, "TableB1_protein_pxd041884.csv")
h_pan, r_pan = read_csv(B1, "TableB1_lactate_axis_panel.csv")
h_ctx, r_ctx = read_csv(B1, "TableB1_context_pdh_pyruvate.csv")
ci = {c: i for i, c in enumerate(h_roll)}
h_roll2 = h_roll + ["lactate_axis_panel_node_class", "note", "in_pdh_pyruvate_context"]
r_roll2 = [r + ["", "", ""] for r in r_roll]
pan_d = as_dict(h_pan, r_pan)
for d in pan_d:
    if d["gene"] in ci:
        pass
gene_to_roll = {r[ci["gene"]]: r for r in r_roll2}
for d in pan_d:
    r = gene_to_roll.get(d["gene"])
    if r is not None:
        r[h_roll2.index("lactate_axis_panel_node_class")] = d["node_class"]
        r[h_roll2.index("note")] = d["note"]
ctx_genes = {row[0] for row in r_ctx}
for r in r_roll2:
    if r[ci["gene"]] in ctx_genes:
        r[h_roll2.index("in_pdh_pyruvate_context")] = "yes"

pidx = {c: i for i, c in enumerate(h_pan)}
roll_uprots = {r[ci["uniprot"]] for r in r_roll2}
app_rows = []
for d in pan_d:
    if d["gene"] in gene_to_roll:
        continue
    row = [d.get(c, "") for c in h_roll]
    row[h_roll.index("test")] = ""
    # a panel entry whose protein is already quantified in the rollup (same UniProt
    # accession) carries the note only; its duplicated quantification cells are left
    # empty because the values live on the rollup entry (e.g., HIST1H4B vs P62805)
    if d["uniprot"] in roll_uprots:
        for c in ["median_DKD", "median_ctrl", "log2FC", "wilcox_p", "percentile_log2FC"]:
            row[h_roll.index(c)] = ""
    app_rows.append(row)
h28, rows28 = concat_by_name([
    (h_roll2, r_roll2, {"block": "protein_rollup_PXD041884"}),
    ([c if c != "node_class" else "lactate_axis_panel_node_class" for c in h_pan], [],
     {"block": "lactate_axis_panel_not_in_rollup"}),
])
# appended rows carry their own panel columns under the union header
pan_ci = {c: i for i, c in enumerate(h_pan)}
for d in pan_d:
    if d["gene"] in gene_to_roll:
        continue
    row = [d.get(h_pan[pan_ci[c]], "") if False else "" for c in h28]
    for c, i in pan_ci.items():
        if c in h28:
            row[h28.index(c)] = d[c]
    row[h28.index("lactate_axis_panel_node_class")] = d["node_class"]
    row[h28.index("test")] = ""
    row[h28.index("in_pdh_pyruvate_context")] = ""
    if d["uniprot"] in roll_uprots:
        for c in ["median_DKD", "median_ctrl", "log2FC", "wilcox_p", "percentile_log2FC"]:
            row[h28.index(c)] = ""
    row[h28.index("block")] = "lactate_axis_panel_not_in_rollup"
    rows28.append(row)
write_csv(os.path.join(BASE, "ESM_28.csv"), h28, rows28)


# ---- ESM_29 / ESM_33: verbatim copies of the module deliverables ----
for src, dst in [(os.path.join(B1, "TableB1_acetylome_inventory.csv"), "ESM_29.csv"),
                 (os.path.join(B3, "TableB3_ilincs_reversers.csv"), "ESM_33.csv")]:
    text = open(src, encoding="utf-8-sig", newline="").read()
    if not text.endswith("\n"):
        text += "\n"
    with open(os.path.join(BASE, dst), "w", encoding="utf-8-sig", newline="") as f:
        f.write(text)


# ---- ESM_30: glomerular iProX panel + global context (explicit column order) ----
h_ip, r_ip = read_csv(B1, "TableB1_iprox_panel.csv")
h_ig, r_ig = read_csv(B1, "TableB1_iprox_global_context.csv")
h30 = ["gene", "detDN_n", "detC_n", "fisher_p", "detection_direction", "medDN", "medC",
       "log2FC_detected_only", "MWU_p", "pct_log2FC", "percentile_log2FC",
       "interpretation_flag", "block"]
rows30 = []
for h, rr, blk in [(h_ip, r_ip, "lactate_axis_panel_glomerular"),
                   (h_ig, r_ig, "global_context_detected_ge50pct_both_groups")]:
    idx = {c: i for i, c in enumerate(h)}
    for r in rr:
        row = [r[idx[c]] if c in idx else "" for c in h30[:-1]] + [blk]
        rows30.append(row)
write_csv(os.path.join(BASE, "ESM_30.csv"), h30, rows30)


# ---- ESM_31: PDK-locus methylation (probes + SUMMARY + region-level + replication) ----
h_t4, r_t4 = read_csv(B4, "TableB4_pdk_methylation.csv")  # run finalize_pdk_summary.py first
h_rl, r_rl = read_csv(B4, "region_level_stats.csv")
h_rp, r_rp = read_csv(B4, "TableB4_replication.csv")
probe_rows = [r for r in r_t4 if r[0].startswith("cg")]
summ_rows = [r for r in r_t4 if not r[0].startswith("cg")]
h31, rows31 = concat_by_name([
    (h_t4, probe_rows, {"block": "probe_GSE233758"}),
    (h_t4, summ_rows, {"block": "summary_GSE233758"}),
    (h_rl, r_rl, {"block": "region_level_GSE233758"}),
    (h_rp, r_rp, {"block": "probe_GSE145745_replication"}),
])
write_csv(os.path.join(BASE, "ESM_31.csv"), h31, rows31)


# ---- ESM_32: LRG-wide methylation sweep, primary + replication ----
h_sw, r_sw = read_csv(B4, "TableB4_lrg_list_methylation.csv")
h_sr, r_sr = read_csv(B4, "TableB4_lrg_replication.csv")
h32, rows32 = concat_by_name([
    (h_sw, r_sw, {"cohort": "GSE233758_primary"}),
    (h_sr, r_sr, {"cohort": "GSE145745_replication"}),
])
write_csv(os.path.join(BASE, "ESM_32.csv"), h32, rows32)

print("ESM_28-33 written to", BASE)
