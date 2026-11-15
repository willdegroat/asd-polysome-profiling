from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
P = HERE.parents[1]
ROOT = P / 'manuscript'
OUT = Path(__file__).resolve().parents[2] / 'manuscript/followup'
d = pd.read_csv(OUT / 'architecture_eligible_analysis_population.csv')
audit = pd.read_csv(P / 'results/transcriptAbundanceComparison.tsv', sep = '\t')
target = 'logThreePrimeUtrLength'
features = [
    'logFivePrimeUtrLength',
    'fivePrimeUtrGc',
    'hasCandidateUorf',
    'hasOverlappingUorf',
    'kozakScore',
    'logCdsLength',
    target,
]
covars = ['logBaseMean', 'logitGeneProportion']
records = []
extra = [
    'inputMeanTpm',
    'inputMeanGeneProportion',
    'inputLog2FoldChange',
    'monosomeLog2FoldChange',
    'polysomeLog2FoldChange',
]
d = d.merge(
    audit[['transcriptId', *extra]], on = 'transcriptId', validate = 'one_to_one'
)
d['logMeanTpm'] = np.log1p(d.meanTpm)
d['logInputMeanTpm'] = np.log1p(d.inputMeanTpm)
p = d.inputMeanGeneProportion.clip(1e-6, 1 - 1e-6)
d['logitInputGeneProportion'] = np.log(p / (1 - p))
d['logGeneProportion'] = np.log(d.meanGeneProportion.clip(1e-6, 1))


def fit(frame, name, predictors = None, outcome = 'log2FoldChange'):
    predictors = features + covars if predictors is None else predictors
    x = frame.dropna(subset = [outcome, *predictors]).copy()
    x = x[x.groupby('parentGeneId').transcriptId.transform('size') >= 2]
    cols = [outcome, *predictors]
    cen = x[cols] - x.groupby('parentGeneId')[cols].transform('mean')
    sd = cen.std(ddof = 1)
    keep = [p for p in predictors if sd[p] > 1e-12]
    X = (cen[keep] / sd[keep]).to_numpy()
    Y = (cen[outcome] / sd[outcome]).to_numpy()
    beta = np.linalg.lstsq(X, Y, rcond = None)[0]
    res = Y - X @ beta
    inv = np.linalg.pinv(X.T @ X)
    ids, codes = np.unique(x.parentGeneId, return_inverse = True)
    scores = np.zeros((len(ids), len(keep)))
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
    j = keep.index(target)
    other = np.delete(X, j, axis = 1)
    xr = X[:, j] - other @ np.linalg.lstsq(other, X[:, j], rcond = None)[0]
    vif = float(np.sum(X[:, j] ** 2) / np.sum(xr**2))
    records.append(
        {
            'model': name,
            'grouping': 'parentGeneId',
            'outcome': outcome,
            'transcripts': n,
            'genes': len(ids),
            'matchingGroups': len(ids),
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
            'targetVIF': vif,
            'designConditionNumber': np.linalg.cond(X),
            'predictors': ';'.join(keep),
        }
    )


for name, cv in [
    ('TPM_and_proportion', ['logMeanTpm', 'logitGeneProportion']),
    ('Input_covariates', ['logInputMeanTpm', 'logitInputGeneProportion']),
    ('log_instead_logit_proportion', ['logBaseMean', 'logGeneProportion']),
]:
    fit(d, 'covariates_' + name, predictors = features + cv)
common = d.dropna(
    subset = [
        'inputLog2FoldChange',
        'monosomeLog2FoldChange',
        'polysomeLog2FoldChange',
        'log2FoldChange',
    ]
)
for layer, col in [
    ('Input', 'inputLog2FoldChange'),
    ('Monosome', 'monosomeLog2FoldChange'),
    ('Polysome', 'polysomeLog2FoldChange'),
    ('PM', 'log2FoldChange'),
]:
    fit(common, 'fitted_common_population_' + layer, outcome = col)
# Six individual-culture profiles on exactly the same well-measured population.
obs = pd.read_csv(OUT / 'observed_PM_fixed_population.csv')
samples = ['07c06', '07c13', '07c16', '1851', '1861', '1916']
ids = set(obs.transcriptId)
q = {}
for fraction, prefix in [
    ('Input', 'input'),
    ('Monosome', 'mono'),
    ('Polysome', 'poly'),
]:
    pieces = []
    for chunk in pd.read_csv(
        P / 'data' / fraction / 'salmon.merged.transcript_tpm.tsv',
        sep = '\t',
        usecols = ['tx'] + [prefix + i for i in samples],
        chunksize = 30000,
    ):
        pieces.append(chunk[chunk.tx.isin(ids)])
    z = pd.concat(pieces).set_index('tx')
    z.columns = samples
    q[fraction] = z
for sample in samples:
    for fraction in ['Monosome', 'Polysome']:
        obs['sampleLogTPM'] = np.log2(
            q[fraction].loc[obs.transcriptId, sample]
        ).to_numpy()
        fit(obs, 'culture_' + sample + '_' + fraction, outcome = 'sampleLogTPM')
    obs['samplePM'] = np.log2(
        q['Polysome'].loc[obs.transcriptId, sample]
        / q['Monosome'].loc[obs.transcriptId, sample]
    ).to_numpy()
    fit(obs, 'culture_' + sample + '_PM', outcome = 'samplePM')
pd.DataFrame(records).to_csv(
    OUT / 'covariates_and_individual_cultures.csv', index = False
)
print(
    pd.DataFrame(records)[
        [
            'model',
            'transcripts',
            'genes',
            'standardizedCoefficient',
            'tStatistic',
            'log2EffectPerDoubledUTR',
        ]
    ].to_string(index = False)
)
