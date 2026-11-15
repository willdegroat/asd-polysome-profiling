import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
P = HERE.parents[1]
ROOT = P / 'manuscript'
OUT = Path(__file__).resolve().parents[2] / 'manuscript/followup'
d = pd.read_csv(OUT / 'architecture_eligible_analysis_population.csv')
seq = json.loads((OUT / 'matched_sequences.json').read_text())['sequences']
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
d['utr3Gc'] = d.transcriptId.map(
    lambda t: (
        (seq[t]['utr3Sequence'].count('G') + seq[t]['utr3Sequence'].count('C'))
        / len(seq[t]['utr3Sequence'])
        if seq[t]['utr3Sequence']
        else np.nan
    )
)
d['sameCDSgroup'] = (
    d.parentGeneId + '|' + d.transcriptId.map(lambda t: seq[t]['cdsHash'])
)
d['sameCDS5group'] = (
    d.sameCDSgroup + '|' + d.transcriptId.map(lambda t: seq[t]['utr5Hash'])
)


def fit(
    frame, name, group = 'parentGeneId', predictors = None, outcome = 'log2FoldChange'
):
    predictors = features + covars if predictors is None else predictors
    x = frame.dropna(subset = [outcome, group, *predictors]).copy()
    x = x[x.groupby(group).transcriptId.transform('size') >= 2]
    if group != 'parentGeneId':
        x = x[x.groupby(group)[target].transform('nunique') >= 2]
    cols = [outcome, *predictors]
    cen = x[cols] - x.groupby(group)[cols].transform('mean')
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
    r = {
        'model': name,
        'grouping': group,
        'outcome': outcome,
        'transcripts': n,
        'genes': len(ids),
        'matchingGroups': x[group].nunique(),
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
        'predictors': ';'.join(keep),
    }
    records.append(r)
    return x


for label, cov in [
    ('none', []),
    ('abundance_only', ['logBaseMean']),
    ('proportion_only', ['logitGeneProportion']),
    ('both', covars),
]:
    fit(d, 'covariates_' + label, predictors = features + cov)
    fit(
        d,
        'sameCDS5_covariates_' + label,
        group = 'sameCDS5group',
        predictors = [target, *cov],
    )
fit(d, 'add_3UTR_GC', predictors = features + covars + ['utr3Gc'])
fit(
    d,
    'sameCDS5_add_3UTR_GC',
    group = 'sameCDS5group',
    predictors = [target, *covars, 'utr3Gc'],
)
for flag in ['gencodePrimary', 'basic']:
    fit(d[d[flag]], 'only_' + flag)
# Exclude groups with very small annotated UTR differences, keeping matched-
# group comparisons.
span = d.groupby('sameCDS5group')[target].agg(lambda a: a.max() - a.min())
ok = set(span[span >= np.log(1.5)].index)
fit(
    d[d.sameCDS5group.isin(ok)],
    'sameCDS5_UTR_ratio_ge_1.5',
    group = 'sameCDS5group',
    predictors = [target, *covars],
)
# Select one large length-contrast pair per gene, without using effects,
# significance or replicate data.
pairs = []
for gid, fr in d.groupby('sameCDS5group', sort = False):
    for a, b in itertools.combinations(list(fr.itertuples()), 2):
        short, long = sorted([a, b], key = lambda a: a.threePrimeUtrLength)
        if (
            short.threePrimeUtrLength < 50
            or long.threePrimeUtrLength - short.threePrimeUtrLength < 100
            or (long.threePrimeUtrLength + 1) / (short.threePrimeUtrLength + 1)
            < 1.5
        ):
            continue
        pairs.append(
            {
                'parentGeneId': long.parentGeneId,
                'geneName': long.geneName,
                'shortId': short.transcriptId,
                'longId': long.transcriptId,
                'shortName': short.transcriptName,
                'longName': long.transcriptName,
                'shortUTR': short.threePrimeUtrLength,
                'longUTR': long.threePrimeUtrLength,
                'logUTRratio': np.log(
                    (long.threePrimeUtrLength + 1)
                    / (short.threePrimeUtrLength + 1)
                ),
                'nested3UTR': seq[long.transcriptId]['utr3Sequence'].startswith(
                    seq[short.transcriptId]['utr3Sequence']
                ),
                'shortPMEffect': short.log2FoldChange,
                'longPMEffect': long.log2FoldChange,
                'longMinusShortPMEffect': long.log2FoldChange
                - short.log2FoldChange,
                'shortAdjustedP': short.adjustedPValue,
                'longAdjustedP': long.adjustedPValue,
                'shortMeanProportion': short.meanGeneProportion,
                'longMeanProportion': long.meanGeneProportion,
                'fuAsd185': long.fuAsd185,
                'synGO': long.synGO,
            }
        )
pairs = (
    pd.DataFrame(pairs)
    .sort_values(
        ['logUTRratio', 'longId', 'shortId'], ascending = [False, True, True]
    )
    .drop_duplicates('parentGeneId')
    .reset_index(drop = True)
)
meta = pd.read_csv(P / 'data/metadata.tsv', sep = '\t', dtype = str)
samples = meta.sampleId.tolist()
ids = set(d.transcriptId)


def read_quant(fraction, prefix, kind = 'tpm'):
    path = P / 'data' / fraction / f'salmon.merged.transcript_{kind}.tsv'
    pieces = []
    for chunk in pd.read_csv(
        path,
        sep = '\t',
        usecols = ['tx'] + [prefix + i for i in samples],
        chunksize = 30000,
    ):
        pieces.append(chunk[chunk.tx.isin(ids)])
    result = pd.concat(pieces).set_index('tx')
    assert set(result.index) == ids
    result.columns = samples
    return result


q = {}
for fraction, prefix in [
    ('Input', 'input'),
    ('Monosome', 'mono'),
    ('Polysome', 'poly'),
]:
    q[fraction] = read_quant(fraction, prefix)
for fraction, prefix in [('Monosome', 'mono'), ('Polysome', 'poly')]:
    q[fraction + 'Counts'] = read_quant(fraction, prefix, 'counts')
# Hold one fixed population for all six omissions; no pseudocounts.
valid = (
    (q['Monosome'].min(axis = 1) >= 1)
    & (q['Polysome'].min(axis = 1) >= 1)
    & (q['MonosomeCounts'].min(axis = 1) >= 10)
    & (q['PolysomeCounts'].min(axis = 1) >= 10)
)
obs = d[d.transcriptId.isin(valid[valid].index)].copy()
logs = np.log2(q['Polysome'].loc[obs.transcriptId]) - np.log2(
    q['Monosome'].loc[obs.transcriptId]
)
obs['observedPM'] = (
    logs.iloc[:, 3:].mean(axis = 1).to_numpy()
    - logs.iloc[:, :3].mean(axis = 1).to_numpy()
)
fit(obs, 'observed_six_culture_PM', outcome = 'observedPM')
for sample in samples:
    left = [s for s in samples[:3] if s != sample]
    right = [s for s in samples[3:] if s != sample]
    obs['omissionPM'] = (
        logs[right].mean(axis = 1).to_numpy() - logs[left].mean(axis = 1).to_numpy()
    )
    fit(obs, 'observed_omit_' + sample, outcome = 'omissionPM')
# Common-population fraction decomposition, with identical predictors and
# transcript rows.
allvalid = valid & (q['Input'].min(axis = 1) >= 1)
common = d[d.transcriptId.isin(allvalid[allvalid].index)].copy()
for fraction in ['Input', 'Monosome', 'Polysome']:
    log_tpm = np.log2(q[fraction].loc[common.transcriptId])
    common['observed' + fraction] = (
        log_tpm.iloc[:, 3:].mean(axis = 1).to_numpy()
        - log_tpm.iloc[:, :3].mean(axis = 1).to_numpy()
    )
    fit(common, 'fraction_component_' + fraction, outcome = 'observed' + fraction)
common['observedPM'] = common.observedPolysome - common.observedMonosome
fit(common, 'fraction_component_PM', outcome = 'observedPM')
# Check examples in the 1:1 same-CDS/same-5UTR selected pairs.
profiles = []
rows = []
for r in pairs.to_dict(orient = 'records'):
    short, long = r['shortId'], r['longId']
    vals = []
    inp = []
    mins = []
    counts = []
    for sample, condition in zip(samples, meta.condition):
        sv = {
            f: float(q[f].loc[short, sample])
            for f in ['Input', 'Monosome', 'Polysome']
        }
        lv = {f: float(q[f].loc[long, sample]) for f in sv}
        z = [sv['Monosome'], lv['Monosome'], sv['Polysome'], lv['Polysome']]
        v = (
            np.log2(lv['Polysome'] / sv['Polysome'])
            - np.log2(lv['Monosome'] / sv['Monosome'])
            if min(z) > 0
            else np.nan
        )
        i = (
            np.log2(lv['Input'] / sv['Input'])
            if min(lv['Input'], sv['Input']) > 0
            else np.nan
        )
        vals.append(v)
        inp.append(i)
        mins.extend(z)
        counts.extend(
            [
                float(q[f + 'Counts'].loc[t, sample])
                for f in ['Monosome', 'Polysome']
                for t in [short, long]
            ]
        )
        profiles.append(
            {
                'geneName': r['geneName'],
                'parentGeneId': r['parentGeneId'],
                'sampleId': sample,
                'condition': condition,
                'longRelativeToShortLog2PM': v,
                'longRelativeToShortLog2Input': i,
                **{f'short{k}TPM': v for k, v in sv.items()},
                **{f'long{k}TPM': v for k, v in lv.items()},
            }
        )
    vals = np.array(vals)
    inp = np.array(inp)
    r.update(
        {
            'minMonosomePolysomeTPM': min(mins),
            'minMonosomePolysomeCount': min(counts),
            'allSixRatiosFinite': bool(np.isfinite(vals).all()),
            'wellMeasuredAllSix': bool(
                np.isfinite(vals).all() and min(mins) >= 1 and min(counts) >= 10
            ),
            'bothMeanProportionGe10pct': min(
                r['shortMeanProportion'], r['longMeanProportion']
            )
            >= 0.1,
        }
    )
    if np.isfinite(vals).all():
        r['observedASDminusUnaffectedLongShortLogPM'] = float(
            vals[3:].mean() - vals[:3].mean()
        )
        loo = [
            np.delete(vals[3:], j - 3).mean() - vals[:3].mean()
            if j >= 3
            else vals[3:].mean() - np.delete(vals[:3], j).mean()
            for j in range(6)
        ]
        r.update(
            {
                'minOmissionPMDifference': min(loo),
                'maxOmissionPMDifference': max(loo),
                'allOmissionsNegative': bool(max(loo) < 0),
                'allOmissionsPositive': bool(min(loo) > 0),
            }
        )
    if np.isfinite(inp).all():
        r['observedASDminusUnaffectedLongShortLogInput'] = float(
            inp[3:].mean() - inp[:3].mean()
        )
    rows.append(r)
pd.DataFrame(rows).to_csv(
    OUT / 'same_CDS_5UTR_large_contrast_pairs.csv', index = False
)
pd.DataFrame(profiles).to_csv(
    OUT / 'same_CDS_5UTR_pair_culture_profiles.csv', index = False
)
pd.DataFrame(records).to_csv(OUT / 'additional_model_results.csv', index = False)
obs.to_csv(OUT / 'observed_PM_fixed_population.csv', index = False)
summ = {
    'largeContrastPairGenes': len(rows),
    'matchedGenes': d[
        d.groupby('sameCDS5group').transcriptId.transform('size') >= 2
    ]
    .groupby('parentGeneId')
    .agg(
        geneName = ('geneName', 'first'),
        fuAsd185 = ('fuAsd185', 'first'),
        synGO = ('synGO', 'first'),
    )
    .reset_index()
    .to_dict(orient = 'records'),
    'observedPopulationCorrelation': float(
        np.corrcoef(obs.log2FoldChange, obs.observedPM)[0, 1]
    ),
}
(OUT / 'additional_summary.json').write_text(json.dumps(summ, indent = 2))
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
print('Large length contrast pairs:', len(rows))
print('Source/observed PM correlation:', summ['observedPopulationCorrelation'])
