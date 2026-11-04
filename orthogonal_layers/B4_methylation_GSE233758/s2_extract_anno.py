"""Extract needed annotation columns from EPIC anno rda files via low-level rdata parser."""
import os
import rdata.parser as rp
import pandas as pd
import numpy as np
import pickle

BASE = BASE
D = BASE + r'\downloads\IlluminaHumanMethylationEPICanno.ilm10b4.hg19\data'


def sym_name(symobj):
    v = symobj.value
    try:
        return v.value.decode()
    except Exception:
        return str(v)


def get_top_object(path):
    r = rp.parse_file(path)
    root = r.object
    # root is a pairlist of (tag=symbol, car=object) chains; take the first car
    return root.value[0] if isinstance(root.value, tuple) else root


def s4_to_dict(obj):
    """Given S4 object, return dict of slot name -> python value/list."""
    slots = {}
    attrs = obj.attributes
    while attrs is not None and getattr(attrs, 'tag', None) is not None:
        name = sym_name(attrs.tag)
        car = attrs.value[0]
        slots[name] = car
        attrs = attrs.value[1] if isinstance(attrs.value, tuple) else None
    return slots


def pairlist_to_dict(obj):
    """Tagged pairlist -> dict of name -> RObject."""
    out = {}
    while obj is not None and getattr(obj, 'tag', None) is not None:
        name = sym_name(obj.tag)
        car = obj.value[0]
        out[name] = car
        obj = obj.value[1] if isinstance(obj.value, tuple) else None
    return out


def rvec_to_list(obj):
    """R STR/CHAR vector -> python list of str (NA -> '')."""
    v = obj.value
    if isinstance(v, list):  # STR vector of CHARs
        res = []
        for ch in v:
            try:
                res.append(ch.value.decode('utf-8', 'replace'))
            except Exception:
                res.append('')
        return res
    if isinstance(v, np.ndarray):
        if v.dtype.kind in 'iuf':
            return v
        return [str(x) for x in v]
    return v


def extract_named_list(obj):
    """VECSXP with names attribute OR tagged pairlist -> dict name -> RObject."""
    out = {}
    t = str(obj.info.type)
    if hasattr(obj.info.type, 'name') and obj.info.type.name == 'VEC' or 'VEC' in t:
        attrs = obj.attributes
        names = None
        while attrs is not None and getattr(attrs, 'tag', None) is not None:
            if sym_name(attrs.tag) == 'names':
                names = rvec_to_list(attrs.value[0])
            attrs = attrs.value[1] if isinstance(attrs.value, tuple) else None
        for i, item in enumerate(obj.value):
            out[names[i] if names else i] = item
        return out
    return pairlist_to_dict(obj)


def extract_file(path, wanted=None):
    top = get_top_object(path)
    slots = s4_to_dict(top)
    rownames = rvec_to_list(slots['rownames'])
    n = len(rownames)
    cols = {}
    coldict = extract_named_list(slots['listData'])
    allnames = list(coldict.keys())
    print(path.split('\\')[-1], 'nrow', n, 'all columns:', allnames)
    for cn in allnames:
        if wanted is not None and cn not in wanted:
            continue
        cols[cn] = rvec_to_list(coldict[cn])
    return pd.DataFrame(cols, index=pd.Index(rownames, name='probe')), allnames


other, other_all = extract_file(D + r'\Other.rda')
loc, loc_all = extract_file(D + r'\Locations.rda')
isl, isl_all = extract_file(D + r'\Islands.UCSC.rda')

anno = other.join(loc, how='left').join(isl, how='left')
anno.to_pickle(BASE + r'\anno_epic_full.pkl')
print('anno shape', anno.shape)
print(anno.head(3).to_string())

# filter PDK2/PDK4
mask = anno['UCSC_RefGene_Name'].fillna('').str.contains('PDK2|PDK4', regex=True)
sub = anno[mask]
print('PDK2/PDK4 probes:', len(sub))
print(sub[['UCSC_RefGene_Name', 'UCSC_RefGene_Group', 'chr', 'pos', 'strand', 'Relation_to_Island']].to_string())
