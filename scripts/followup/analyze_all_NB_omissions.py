from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = Path(__file__).resolve().parents[2] / 'manuscript/followup'
d = pd.read_csv(OUT / 'architecture_eligible_analysis_population.csv')
features = [
    'logFivePrimeUtrLength',
    'fivePrimeUtrGc',
    'hasCandidateUorf',
    'hasOverlappingUorf',
    'kozakScore',
    'logCdsLength',
    'logThreePrimeUtrLength',
]
records = []
omits = ['none', '07c06', '07c13', '07c16', '1851', '1861', '1916']
fits = {
    s: pd.read_csv(OUT / f'fixed_dispersion_NB_omit_{s}.csv')
    .set_index('transcriptId')
    .reindex(d.transcriptId)
    for s in omits
}
valid = np.ones(len(d), dtype = bool)
for fit in fits.values():
    valid &= (
        fit.betaConverged & fit.log2FoldChange.notna() & fit.stat.notna()
    ).to_numpy()
x = d.loc[valid].copy()
well = set(pd.read_csv(OUT / 'observed_PM_fixed_population.csv').transcriptId)


def architecture(frame, label, outcome):
    cols = [outcome, *features, 'logBaseMean', 'logitGeneProportion']
    z = frame.dropna(subset = cols).copy()
    z = z[z.groupby('parentGeneId').transcriptId.transform('size') >= 2]
    cen = z[cols] - z.groupby('parentGeneId')[cols].transform('mean')
    sd = cen.std(ddof = 1)
    X = (cen[cols[1:]] / sd[cols[1:]]).to_numpy()
    Y = (cen[outcome] / sd[outcome]).to_numpy()
    beta = np.linalg.lstsq(X, Y, rcond = None)[0]
    res = Y - X @ beta
    inv = np.linalg.pinv(X.T @ X)
    ids, codes = np.unique(z.parentGeneId, return_inverse = True)
    scores = np.zeros((len(ids), X.shape[1]))
    np.add.at(scores, codes, X * res[:, None])
    n, k = X.shape
    V = (
        inv
        @ (scores.T @ scores)
        @ inv
        * len(ids)
        / (len(ids) - 1)
        * (n - 1)
        / (n - k)
    )
    se = np.sqrt(np.maximum(np.diag(V), 0))
    for target in ['logThreePrimeUtrLength', 'logCdsLength']:
        j = cols[1:].index(target)
        records.append(
            {
                'model': label + '_' + target,
                'outcome': outcome,
                'feature': target,
                'transcripts': n,
                'genes': len(ids),
                'standardizedCoefficient': beta[j],
                'standardError': se[j],
                'tStatistic': beta[j] / se[j],
                'df': len(ids) - 1,
                'outcomeUnitsPerSD': beta[j] * sd[outcome],
                'outcomeSD': sd[outcome],
                'featureSD': sd[target],
                'predictors': ';'.join(cols[1:]),
            }
        )


for pop, fr in [
    ('full_frozen', x),
    ('well_measured_frozen', x[x.transcriptId.isin(well)]),
]:
    for omit in omits:
        for kind, col in [('FC', 'log2FoldChange'), ('Wald', 'stat')]:
            z = fr.copy()
            z['NB' + kind] = fits[omit].loc[z.transcriptId, col].to_numpy()
            architecture(z, pop + '_omit_' + omit + '_NB' + kind, 'NB' + kind)
r = pd.DataFrame(records)
r.to_csv(OUT / 'all_omissions_NB_architecture_models.csv', index = False)
pd.DataFrame(
    {'transcriptId': d.transcriptId, 'converged_in_all_seven_fits': valid}
).to_csv(OUT / 'NB_frozen_population_flags.csv', index = False)
print(
    r[r.outcome.eq('NBWald')][
        [
            'model',
            'transcripts',
            'genes',
            'standardizedCoefficient',
            'tStatistic',
        ]
    ].to_string(index = False)
)
