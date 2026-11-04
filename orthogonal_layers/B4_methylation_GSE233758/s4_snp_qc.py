"""Flag common SNPs at the 61 candidate probes using SNPs.147CommonSingle.rda."""
import os
import rdata.parser as rp
import pandas as pd
import numpy as np

BASE = BASE
D = BASE + r'\downloads\IlluminaHumanMethylationEPICanno.ilm10b4.hg19\data'
p = D + r'\SNPs.147CommonSingle.rda'

def sym_name(symobj):
    v = symobj.value
    try:
        return v.value.decode()
    except Exception:
        return str(v)

def get_top_object(path):
    r = rp.parse_file(path)
    root = r.object
    return root.value[0] if isinstance(root.value, tuple) else root

def s4_to_dict(obj):
    slots = {}
    attrs = obj.attributes
    while attrs is not None and getattr(attrs, 'tag', None) is not None:
        slots[sym_name(attrs.tag)] = attrs.value[0]
        attrs = attrs.value[1] if isinstance(attrs.value, tuple) else None
    return slots

def rvec_to_list(obj):
    v = obj.value
    if isinstance(v, list):
        out = []
        for ch in v:
            cv = getattr(ch, 'value', None)
            if cv is None:
                out.append('')
            elif isinstance(cv, bytes):
                out.append(cv.decode('utf-8', 'replace'))
            else:
                out.append(str(cv))
        return out
    if isinstance(v, np.ndarray):
        return v
    return v

def extract_named_list(obj):
    out = {}
    t = str(obj.info.type)
    if 'VEC' in t:
        attrs = obj.attributes
        names = None
        while attrs is not None and getattr(attrs, 'tag', None) is not None:
            if sym_name(attrs.tag) == 'names':
                names = rvec_to_list(attrs.value[0])
            attrs = attrs.value[1] if isinstance(attrs.value, tuple) else None
        for i, item in enumerate(obj.value):
            out[names[i] if names else i] = item
        return out
    # pairlist
    while obj is not None and getattr(obj, 'tag', None) is not None:
        out[sym_name(obj.tag)] = obj.value[0]
        obj = obj.value[1] if isinstance(obj.value, tuple) else None
    return out

top = get_top_object(p)
slots = s4_to_dict(top)
rownames = rvec_to_list(slots['rownames'])
coldict = extract_named_list(slots['listData'])
print('SNP columns:', list(coldict.keys()))
cols = {}
for cn, co in coldict.items():
    cols[cn] = rvec_to_list(co)
snp = pd.DataFrame(cols, index=pd.Index(rownames, name='probe'))

cand = pd.read_csv(BASE + r'\probe_stats_raw.csv')
flagged = []
for pid in cand['cgID']:
    if pid in snp.index:
        r = snp.loc[pid]
        hits = [f'{c}={r[c]}' for c in snp.columns if str(r[c]) not in ('', 'nan', 'None')]
        if any('' != str(r[c]) and str(r[c]) != 'nan' for c in ['Probe_rs', 'CpG_rs', 'SBE_rs']):
            flagged.append((pid, '; '.join(hits)))
print('probes with common SNP at CpG/SBE/probe:')
for f in flagged:
    print(' ', f)
pd.DataFrame(flagged, columns=['probe', 'snp_hits']).to_csv(BASE + r'\snp_flags.csv', index=False)
