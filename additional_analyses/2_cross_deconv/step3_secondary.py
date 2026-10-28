import pandas as pd
from scipy.stats import mannwhitneyu
F = pd.read_csv('evidence_gain/2_cross_deconv/deconv_fractions_all.csv')
c = F[F.cohort == 'GSE142025']
ea = c[c.group.isin(['Advanced_DN', 'Early_DN'])]
adv, ear = ea[ea.group == 'Advanced_DN'], ea[ea.group == 'Early_DN']
u = mannwhitneyu(adv['B'], ear['B'], alternative='two-sided')
print(f"GSE142025 B fraction Adv(n={len(adv)}) vs Early(n={len(ear)}): median {adv['B'].median():.3f} vs {ear['B'].median():.3f}, p={u.pvalue:.4f}")
co = pd.read_csv('evidence_gain/2_cross_deconv/deconv_pdk_correlations.csv')
for feat in ('signature_score', 'PDK4'):
    ss = co[(co.feature == feat) & (co.type == 'B') & (co.cohort.isin(['GSE96804', 'GSE142025']))]
    print(f'\n{feat} vs B fraction:')
    print(ss[['cohort', 'rho_all', 'p_all', 'rho_within_DN']].to_string(index=False))
