import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
P = HERE.parents[1]
ROOT = P / 'manuscript'
OUT = Path(__file__).resolve().parents[2] / 'manuscript/followup'
OUT.mkdir(exist_ok = True)

features = [
    'logFivePrimeUtrLength',
    'fivePrimeUtrGc',
    'hasCandidateUorf',
    'hasOverlappingUorf',
    'kozakScore',
    'logCdsLength',
    'logThreePrimeUtrLength',
]
covars = ['logBaseMean', 'logitGeneProportion']
target = 'logThreePrimeUtrLength'
files = [
    'transcriptRedistribution.tsv',
    'transcriptFeatures.tsv',
    'transcriptAbundanceComparison.tsv',
    'transcriptFeatureSensitivity.tsv',
    'fuAsdGenes.tsv',
    'polysomeMonosome.tsv',
]
hashes = {
    n: hashlib.sha256((P / 'results' / n).read_bytes()).hexdigest()
    for n in files
}
tx = pd.read_csv(P / 'results' / files[0], sep = '\t')
f = pd.read_csv(P / 'results' / files[1], sep = '\t')
audit = pd.read_csv(P / 'results' / files[2], sep = '\t')
sens = pd.read_csv(P / 'results' / files[3], sep = '\t')
risk = pd.read_csv(P / 'results' / files[4], sep = '\t')
gene = pd.read_csv(P / 'results' / files[5], sep = '\t')
d = tx.merge(
    f.drop(columns = ['parentGeneId', 'geneName']),
    on = 'transcriptId',
    validate = 'one_to_one',
).merge(
    audit[
        [
            'transcriptId',
            'redistributionStructuralZero',
            'redistributionExtremeEffect',
        ]
    ],
    on = 'transcriptId',
    how = 'left',
    validate = 'one_to_one',
)
d = d[
    d.architectureEligible
    & ~d.redistributionStructuralZero.fillna(False).astype(bool)
    & ~d.redistributionExtremeEffect.fillna(False).astype(bool)
].copy()
d['logBaseMean'] = np.log1p(d.baseMean)
pr = d.meanGeneProportion.clip(1e-6, 1 - 1e-6)
d['logitGeneProportion'] = np.log(pr / (1 - pr))
d = d.dropna(subset = ['parentGeneId', 'log2FoldChange', *features, *covars])
d = d[d.groupby('parentGeneId').transcriptId.transform('size') >= 2].copy()
assert len(d) == 11679 and d.parentGeneId.nunique() == 3843
riskset = set(risk.loc[risk.fuAsd185, 'geneId'])
d['fuAsd185'] = d.parentGeneId.str.split('.').str[0].isin(riskset)
synsets = {}
for line in (P / 'data/geneSets/synGO2024.gmt').read_text().splitlines():
    a = line.split('\t')
    synsets[a[0]] = set(a[2:])
synunion = set().union(*synsets.values())
d['synGO'] = d.geneName.isin(synunion)
records = []


def fit(frame, name, group = 'parentGeneId', predictors = None, interaction = None):
    predictors = features + covars if predictors is None else predictors
    x = frame.dropna(subset = ['log2FoldChange', group, *predictors]).copy()
    x = x[x.groupby(group).transcriptId.transform('size') >= 2]
    if group != 'parentGeneId':
        x = x[x.groupby(group)[target].transform('nunique') >= 2]
    cols = ['log2FoldChange', *predictors]
    cen = x[cols] - x.groupby(group)[cols].transform('mean')
    sd = cen.std(ddof = 1)
    keep = [p for p in predictors if sd[p] > 1e-12]
    Z = cen[keep] / sd[keep]
    Y = (cen.log2FoldChange / sd.log2FoldChange).to_numpy()
    X = Z.to_numpy()
    names = keep.copy()
    if interaction is not None:
        X = np.column_stack(
            [X, Z[target].to_numpy() * x[interaction].astype(float).to_numpy()]
        )
        names.append('UTR3_x_' + interaction)
    beta = np.linalg.lstsq(X, Y, rcond = None)[0]
    res = Y - X @ beta
    inv = np.linalg.pinv(X.T @ X)
    ids, codes = np.unique(x.parentGeneId, return_inverse = True)
    scores = np.zeros((len(ids), len(names)))
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
    j = names.index('UTR3_x_' + interaction if interaction else target)
    r = {
        'model': name,
        'grouping': group,
        'transcripts': n,
        'genes': len(ids),
        'matchingGroups': x[group].nunique(),
        'standardizedCoefficient': beta[j],
        'standardError': se[j],
        'tStatistic': beta[j] / se[j],
        'df': len(ids) - 1,
        'outcomeUnitsPerSD': beta[j] * sd.log2FoldChange,
        'outcomeSD': sd.log2FoldChange,
        'featureSD': sd[target],
        'predictors': ';'.join(names),
    }
    if interaction:
        main = beta[names.index(target)]
        r.update(
            {
                'modifier': interaction,
                'modifierGenes': x.loc[
                    x[interaction], 'parentGeneId'
                ].nunique(),
                'otherSlopeStandardized': main,
                'modifierSlopeStandardized': main + beta[j],
            }
        )
    records.append(r)
    return x, X, Y, beta, res, names, sd


x, X, Y, beta, res, names, sd = fit(d, 'baseline_fold_change')
stored = sens[
    sens.analysis.eq('redistribution')
    & sens.variant.eq('effectSize')
    & sens.model.eq('joint')
].set_index('feature')
assert np.allclose(
    beta[:7], stored.loc[features].standardizedCoefficient, atol = 1e-10
)
for threshold in [0.05, 0.10, 0.20]:
    fit(d[d.meanGeneProportion >= threshold], f'mean_proportion_ge_{threshold}')
for threshold in [50, 100, 500]:
    fit(d[d.baseMean >= threshold], f'baseMean_ge_{threshold}')
fit(d, 'gene_balanced_placeholder') if False else None
# Calculate gene contributions to the adjusted covariance for 3-prime UTR.
j = names.index(target)
other = np.delete(X, j, axis = 1)
xr = X[:, j] - other @ np.linalg.lstsq(other, X[:, j], rcond = None)[0]
yr = Y - other @ np.linalg.lstsq(other, Y, rcond = None)[0]
con = (
    pd.DataFrame(
        {
            'parentGeneId': x.parentGeneId,
            'geneName': x.geneName,
            'covariance': xr * yr,
        }
    )
    .groupby(['parentGeneId', 'geneName'])
    .covariance.sum()
    .reset_index()
)
con['absoluteCovariance'] = con.covariance.abs()
con = con.sort_values('absoluteCovariance', ascending = False)
con['fractionAbsoluteCovariance'] = (
    con.absoluteCovariance / con.absoluteCovariance.sum()
)
con.to_csv(OUT / 'gene_contributions.csv', index = False)
for top in [1, 10, 50]:
    fit(
        d[~d.parentGeneId.isin(con.head(top).parentGeneId)],
        f'exclude_top_{top}_contributing_genes',
    )
for flag in ['fuAsd185', 'synGO']:
    fit(d, 'interaction_' + flag, interaction = flag)
    fit(d[~d[flag]], 'exclude_' + flag)
print(
    'Baseline reproduced; sensitivity/interaction models computed.', flush = True
)
# Exact sequence matching from the run-matched transcript FASTA, cached under
# this task only.
seqcache = OUT / 'matched_sequences.json'
fasta = P / 'data/reference/genome.transcripts.fa'
fst = fasta.stat()
identity = {
    'bytes': fst.st_size,
    'mtime': fst.st_mtime_ns,
    'idsSHA': hashlib.sha256(
        '\n'.join(sorted(d.transcriptId)).encode()
    ).hexdigest(),
}
if (
    seqcache.exists()
    and json.loads(seqcache.read_text()).get('identity') == identity
):
    seq = json.loads(seqcache.read_text())['sequences']
else:
    seq = {}
    metadata = d.set_index('transcriptId')
    targets = set(metadata.index)
    current = None
    parts = []

    def save():
        if current in targets:
            s = ''.join(parts).upper()
            row = metadata.loc[current]
            start = int(row.mainStartIndex)
            stop = int(row.mainStopIndexExclusive)
            assert len(s) == int(row.sequenceLength) and len(s) - stop == int(
                row.threePrimeUtrLength
            )
            seq[current] = {
                'cdsHash': hashlib.sha256(s[start:stop].encode()).hexdigest(),
                'utr5Hash': hashlib.sha256(s[:start].encode()).hexdigest(),
                'utr3Sequence': s[stop:],
            }

    with fasta.open(encoding = 'ascii') as stream:
        for line in stream:
            if line.startswith('>'):
                save()
                current = line[1:].strip().split()[0]
                parts = []
            elif current in targets:
                parts.append(line.strip())
        save()
    assert set(seq) == targets
    seqcache.write_text(json.dumps({'identity': identity, 'sequences': seq}))
print(
    'Exact CDS/5-prime UTR sequences verified for',
    len(seq),
    'transcripts.',
    flush = True,
)
d['cdsHash'] = d.transcriptId.map(lambda t: seq[t]['cdsHash'])
d['utr5Hash'] = d.transcriptId.map(lambda t: seq[t]['utr5Hash'])
d['sameCDSgroup'] = d.parentGeneId + '|' + d.cdsHash
d['sameCDS5group'] = d.sameCDSgroup + '|' + d.utr5Hash
fit(
    d,
    'same_CDS_sequence',
    group = 'sameCDSgroup',
    predictors = [
        target,
        'logFivePrimeUtrLength',
        'fivePrimeUtrGc',
        'hasCandidateUorf',
        'hasOverlappingUorf',
        'kozakScore',
        *covars,
    ],
)
fit(
    d,
    'same_CDS_and_5UTR_sequence',
    group = 'sameCDS5group',
    predictors = [target, *covars],
)
# One nested long/short 3-prime UTR pair per gene, selected without outcome information.
pairs = []
for gid, frame in d.groupby('sameCDS5group', sort = False):
    rows = list(frame.itertuples())
    for a, b in itertools.combinations(rows, 2):
        short, long = sorted([a, b], key = lambda r: r.threePrimeUtrLength)
        if (
            short.threePrimeUtrLength < 50
            or long.threePrimeUtrLength - short.threePrimeUtrLength < 500
            or (long.threePrimeUtrLength + 1) / (short.threePrimeUtrLength + 1)
            < 1.5
        ):
            continue
        if not seq[long.transcriptId]['utr3Sequence'].startswith(
            seq[short.transcriptId]['utr3Sequence']
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
pairs.to_csv(OUT / 'nested_same_CDS_5UTR_pairs.csv', index = False)
print(
    'Nested exact-CDS/5UTR pairs selected for', len(pairs), 'genes.', flush = True
)
d.drop(columns = ['cdsHash', 'utr5Hash', 'sameCDSgroup', 'sameCDS5group']).to_csv(
    OUT / 'architecture_eligible_analysis_population.csv', index = False
)
pd.DataFrame(records).to_csv(OUT / 'new_model_results.csv', index = False)
# Read full quantification tables in chunks; retain only matched-pair transcripts.
meta = pd.read_csv(P / 'data/metadata.tsv', sep = '\t', dtype = str)
samples = meta.sampleId.tolist()
ids = set(pairs.shortId) | set(pairs.longId)


def read_quant(fraction, prefix, kind):
    path = P / 'data' / fraction / f'salmon.merged.transcript_{kind}.tsv'
    columns = ['tx'] + [prefix + i for i in samples]
    pieces = []
    for chunk in pd.read_csv(path, sep = '\t', usecols = columns, chunksize = 30000):
        pieces.append(chunk[chunk.tx.isin(ids)])
    result = pd.concat(pieces).set_index('tx')
    assert set(result.index) == ids
    result.columns = samples
    return result


quant = {}
for fraction, prefix in [
    ('Input', 'input'),
    ('Monosome', 'mono'),
    ('Polysome', 'poly'),
]:
    quant[fraction] = read_quant(fraction, prefix, 'tpm')
    if fraction != 'Input':
        quant[fraction + 'Counts'] = read_quant(fraction, prefix, 'counts')
print('Six-culture matched-pair quantification extracted.', flush = True)
rows = []
profile = []
for r in pairs.to_dict(orient = 'records'):
    short, long = r['shortId'], r['longId']
    vals = []
    inputvals = []
    minimum_tpm = []
    minimum_counts = []
    for sample, condition in zip(samples, meta.condition):
        sl = {
            f: float(quant[f].loc[short, sample])
            for f in ['Input', 'Monosome', 'Polysome']
        }
        lo = {f: float(quant[f].loc[long, sample]) for f in sl}
        pos = all(
            v > 0
            for v in [
                sl['Monosome'],
                lo['Monosome'],
                sl['Polysome'],
                lo['Polysome'],
            ]
        )
        v = (
            np.log2(lo['Polysome'] / sl['Polysome'])
            - np.log2(lo['Monosome'] / sl['Monosome'])
            if pos
            else np.nan
        )
        inp = (
            np.log2(lo['Input'] / sl['Input'])
            if lo['Input'] > 0 and sl['Input'] > 0
            else np.nan
        )
        vals.append(v)
        inputvals.append(inp)
        minimum_tpm.extend(
            [sl['Monosome'], lo['Monosome'], sl['Polysome'], lo['Polysome']]
        )
        minimum_counts.extend(
            [
                float(quant[f + 'Counts'].loc[t, sample])
                for f in ['Monosome', 'Polysome']
                for t in [short, long]
            ]
        )
        profile.append(
            {
                'geneName': r['geneName'],
                'parentGeneId': r['parentGeneId'],
                'sampleId': sample,
                'condition': condition,
                'longRelativeToShortLog2PM': v,
                'longRelativeToShortLog2Input': inp,
                **{f'short{k}TPM': v for k, v in sl.items()},
                **{f'long{k}TPM': v for k, v in lo.items()},
            }
        )
    vals = np.array(vals)
    inp = np.array(inputvals)
    r['minMonosomePolysomeTPM'] = min(minimum_tpm)
    r['minMonosomePolysomeCount'] = min(minimum_counts)
    r['allSixRatiosFinite'] = bool(np.isfinite(vals).all())
    r['wellMeasuredAllSix'] = bool(
        r['allSixRatiosFinite']
        and min(minimum_tpm) >= 1
        and min(minimum_counts) >= 10
    )
    r['bothMeanProportionGe10pct'] = (
        min(r['shortMeanProportion'], r['longMeanProportion']) >= 0.1
    )
    if np.isfinite(vals).all():
        r['UnaffectedMeanLongShortLogPM'] = float(vals[:3].mean())
        r['ASDMeanLongShortLogPM'] = float(vals[3:].mean())
        r['observedASDminusUnaffectedLongShortLogPM'] = float(
            vals[3:].mean() - vals[:3].mean()
        )
        loo = [
            np.delete(vals[3:], j - 3).mean() - vals[:3].mean()
            if j >= 3
            else vals[3:].mean() - np.delete(vals[:3], j).mean()
            for j in range(6)
        ]
        r['minOmissionPMDifference'] = float(min(loo))
        r['maxOmissionPMDifference'] = float(max(loo))
        r['allOmissionsNegative'] = bool(max(loo) < 0)
        r['allOmissionsPositive'] = bool(min(loo) > 0)
    if np.isfinite(inp).all():
        r['observedASDminusUnaffectedLongShortLogInput'] = float(
            inp[3:].mean() - inp[:3].mean()
        )
    rows.append(r)
pairs = pd.DataFrame(rows)
profiles = pd.DataFrame(profile)
pairs.to_csv(OUT / 'nested_pairs_with_replicate_diagnostics.csv', index = False)
profiles.to_csv(OUT / 'nested_pairs_six_culture_profiles.csv', index = False)
summary = {
    'baselineGenes': int(d.parentGeneId.nunique()),
    'baselineTranscripts': len(d),
    'ASD185genes': int(d.loc[d.fuAsd185, 'parentGeneId'].nunique()),
    'SynGOgenes': int(d.loc[d.synGO, 'parentGeneId'].nunique()),
    'topGeneCovarianceFraction': float(con.fractionAbsoluteCovariance.iloc[0]),
    'nestedPairGenes': len(pairs),
}
for name, frame in [
    ('all_finite', pairs[pairs.allSixRatiosFinite]),
    ('well_measured', pairs[pairs.wellMeasuredAllSix]),
    (
        'well_measured_major_isoforms',
        pairs[pairs.wellMeasuredAllSix & pairs.bothMeanProportionGe10pct],
    ),
]:
    if len(frame):
        values = frame.observedASDminusUnaffectedLongShortLogPM
        omissions = {}
        for sample in samples:
            pr = profiles[
                profiles.parentGeneId.isin(frame.parentGeneId)
                & profiles.sampleId.ne(sample)
            ]
            avg = (
                pr.groupby(['parentGeneId', 'condition'])
                .longRelativeToShortLog2PM.mean()
                .unstack()
            )
            omissions[sample] = float((avg.ASD - avg.Unaffected).median())
        summary[name] = {
            'genes': len(frame),
            'negativeGenes': int((values < 0).sum()),
            'fractionNegative': float((values < 0).mean()),
            'medianLog2Difference': float(values.median()),
            'medianFoldRatio': float(np.exp2(values.median())),
            'omissionMedianDifferences': omissions,
            'allOmissionsNegativeGenes': int(frame.allOmissionsNegative.sum()),
            'ASD185genes': int(frame.fuAsd185.sum()),
            'SynGOgenes': int(frame.synGO.sum()),
        }
summary['sequenceFASTAIdentity'] = identity
summary['sourceSHA256'] = hashes
(OUT / 'summary.json').write_text(json.dumps(summary, indent = 2))
for n, h in hashes.items():
    assert hashlib.sha256((P / 'results' / n).read_bytes()).hexdigest() == h
print(json.dumps(summary, indent = 2))
