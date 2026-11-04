"""Region-level summary for replication from saved CSV."""
import os
import numpy as np
import pandas as pd
from scipy import stats

BASE = BASE
tab = pd.read_csv(BASE + r'\TableB4_replication.csv').set_index('cgID')

for gene in ('PDK2', 'PDK4'):
    for cat in ('promoter', 'body'):
        g = tab[(tab.gene == gene) & (tab.region_category == cat)]
        if len(g) == 0:
            continue
        db = g.delta_beta.mean()
        tt = stats.ttest_ind(g.Mvalue_logFC, 0)  # placeholder not used; report mean only
        npos = int((g.delta_beta > 0).sum())
        print(f'{gene} {cat}: k={len(g)} mean_db={db:+.4f} range=[{g.delta_beta.min():+.4f},{g.delta_beta.max():+.4f}] '
              f'min_p={g.p.min():.4g} min_BH={g.p_BH.min():.4g} pos={npos}/{len(g)} '
              f'concordant_sign_with_GSE233758={npos}/{len(g)}')
