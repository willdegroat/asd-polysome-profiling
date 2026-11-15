import json
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
a = pd.read_csv(OUT / 'fixed_dispersion_NB_omit_none.csv').set_index(
    'transcriptId'
)
b = pd.read_csv(OUT / 'fixed_dispersion_NB_omit_07c16.csv').set_index(
    'transcriptId'
)
a = a.reindex(d.transcriptId).reset_index()
b = b.reindex(d.transcriptId).reset_index()
comparison = {
    'allTranscripts': len(d),
    'fullFitNonconverged': int((~a.betaConverged).sum()),
    'omitFitNonconverged': int((~b.betaConverged).sum()),
    'fullVsSourceLFCcorrelation': float(
        np.corrcoef(d.log2FoldChange, a.log2FoldChange)[0, 1]
    ),
    'medianAbsoluteFullVsSourceLFCdifference': float(
        np.median(abs(d.log2FoldChange - a.log2FoldChange))
    ),
    'p95AbsoluteFullVsSourceLFCdifference': float(
        np.quantile(abs(d.log2FoldChange - a.log2FoldChange), 0.95)
    ),
}
valid = (
    a.betaConverged
    & b.betaConverged
    & a.log2FoldChange.notna()
    & b.log2FoldChange.notna()
    & a.stat.notna()
    & b.stat.notna()
).to_numpy()
x = d.loc[valid].copy()
x['sourceFC'] = x.log2FoldChange
x['sourceWald'] = x.statistic
for key, col in [('FC', 'log2FoldChange'), ('Wald', 'stat')]:
    x['fullNB' + key] = a.loc[valid, col].to_numpy()
    x['omit07c16NB' + key] = b.loc[valid, col].to_numpy()
obs = pd.read_csv(OUT / 'observed_PM_fixed_population.csv')
well = set(obs.transcriptId)


def fit(frame, name, outcome):
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
                'model': name + '_' + target,
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
                'log2EffectPerDoubledUTR': beta[j]
                * sd[outcome]
                / sd[target]
                * np.log(2),
                'predictors': ';'.join(cols[1:]),
            }
        )


for pop, fr in [
    ('full_frozen', x),
    ('well_measured_frozen', x[x.transcriptId.isin(well)]),
]:
    for outcome in [
        'sourceFC',
        'sourceWald',
        'fullNBFC',
        'fullNBWald',
        'omit07c16NBFC',
        'omit07c16NBWald',
    ]:
        fit(fr, pop + '_' + outcome, outcome)
pd.DataFrame(records).to_csv(
    OUT / 'fixed_dispersion_NB_architecture_models.csv', index = False
)
(OUT / 'NB_refit_reproduction.json').write_text(
    json.dumps(comparison, indent = 2)
)
print(json.dumps(comparison, indent = 2))
print(
    pd.DataFrame(records)[
        [
            'model',
            'transcripts',
            'genes',
            'standardizedCoefficient',
            'tStatistic',
        ]
    ].to_string(index = False)
)
