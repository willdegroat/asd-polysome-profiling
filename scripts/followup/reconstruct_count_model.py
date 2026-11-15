from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
P = HERE.parents[1]
OUT = Path(__file__).resolve().parents[2] / 'manuscript/followup'
tx = pd.read_csv(P / 'results/transcriptRedistribution.tsv', sep = '\t')
ids = set(tx.transcriptId)
samples = ['07c06', '07c13', '07c16', '1851', '1861', '1916']
allcounts = []
for fraction, folder, prefix in [
    ('input', 'Input', 'input'),
    ('monosome', 'Monosome', 'mono'),
    ('polysome', 'Polysome', 'poly'),
]:
    names = ['tx'] + [prefix + s for s in samples]
    streams = [
        pd.read_csv(
            P / 'data' / folder / f'salmon.merged.transcript_{kind}.tsv',
            sep = '\t',
            usecols = names,
            chunksize = 30000,
        )
        for kind in ['counts', 'tpm', 'lengths']
    ]
    cs = np.zeros(6)
    ss = np.zeros(6)
    pieces = []
    for c, t, lengths in zip(*streams):
        assert c.tx.equals(t.tx) and c.tx.equals(lengths.tx)
        scaled = (
            t.iloc[:, 1:].to_numpy()
            * lengths.iloc[:, 1:].mean(axis = 1).to_numpy()[:, None]
        )
        cs += c.iloc[:, 1:].sum(axis = 0).to_numpy()
        ss += scaled.sum(axis = 0)
        keep = c.tx.isin(ids)
        z = pd.DataFrame(scaled[keep], index = c.loc[keep, 'tx'], columns = samples)
        pieces.append(z)
    z = pd.concat(pieces).mul(cs / ss, axis = 1).reindex(tx.transcriptId)
    assert z.notna().all().all()
    z = z.round().astype('int64')
    z.columns = [s + '_' + fraction for s in samples]
    allcounts.append(z)
    print(
        'Reconstructed',
        fraction,
        'counts for',
        len(z),
        'transcripts.',
        flush = True,
    )
counts = pd.concat(allcounts, axis = 1)
counts.index.name = 'transcriptId'
counts.to_csv(OUT / 'count_model_reconstructed_counts.csv')
print('Ready for fixed-dispersion NB sensitivity refits.', flush = True)
