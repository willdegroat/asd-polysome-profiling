"""Reproduce the exact-sequence-matched 3â€² UTR analysis in Figure 3D."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[1]
TARGET = 'logThreePrimeUtrLength'


def fit_model(
    table: pd.DataFrame, *, adjusted: bool
) -> tuple[dict, pd.DataFrame]:
    """Fit within-matching-group effects with parent-gene-clustered uncertainty."""
    predictors = [TARGET]
    if adjusted:
        predictors += ['logBaseMean', 'logitGeneProportion']
    variables = ['log2FoldChange', *predictors]
    centered = table[variables] - table.groupby('matchingGroup')[
        variables
    ].transform('mean')
    sd = centered.std(ddof = 1)
    x = (centered[predictors] / sd[predictors]).to_numpy()
    y = (centered.log2FoldChange / sd.log2FoldChange).to_numpy()
    beta = np.linalg.lstsq(x, y, rcond = None)[0]
    residual = y - x @ beta
    genes, codes = np.unique(table.parentGeneId, return_inverse = True)
    scores = np.zeros((len(genes), len(predictors)))
    np.add.at(scores, codes, x * residual[:, None])
    inverse = np.linalg.pinv(x.T @ x)
    n, k = x.shape
    correction = len(genes) / (len(genes) - 1) * (n - 1) / (n - k)
    covariance = inverse @ scores.T @ scores @ inverse * correction
    se = np.sqrt(covariance[0, 0])
    df = len(genes) - 1
    critical = t.ppf(0.975, df)
    record = {
        'model': 'adjusted' if adjusted else 'unadjusted',
        'transcripts': n,
        'genes': len(genes),
        'matchingGroups': table.matchingGroup.nunique(),
        'standardizedCoefficient': beta[0],
        'standardError': se,
        'pValueUnadjusted': 2 * t.sf(abs(beta[0] / se), df),
        'confidenceIntervalLow': beta[0] - critical * se,
        'confidenceIntervalHigh': beta[0] + critical * se,
    }
    points = table.copy()
    points['withinMatchedGroupUTRlengthSD'] = x[:, 0]
    points['adjustedPMpartialResidualLog2'] = (
        residual + beta[0] * x[:, 0]
    ) * sd.log2FoldChange
    return record, points


def main() -> None:
    """Reconstruct matching groups and export adjusted and unadjusted models."""
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument(
        '--output', type = Path, default = ROOT / 'manuscript/work/matched_utr3'
    )
    args = parser.parse_args()
    folder = ROOT / 'manuscript/followup'
    table = pd.read_csv(
        folder / 'architecture_eligible_analysis_population.csv'
    )
    sequences = json.loads((folder / 'matched_sequences.json').read_text())[
        'sequences'
    ]
    table['matchingGroup'] = (
        table.parentGeneId
        + '|'
        + table.transcriptId.map(lambda key: sequences[key]['cdsHash'])
        + '|'
        + table.transcriptId.map(lambda key: sequences[key]['utr5Hash'])
    )
    table = table[
        table.groupby('matchingGroup').transcriptId.transform('size') >= 2
    ]
    table = table[
        table.groupby('matchingGroup')[TARGET].transform('nunique') >= 2
    ]
    assert (len(table), table.parentGeneId.nunique()) == (154, 69)
    records = []
    args.output.mkdir(parents = True, exist_ok = True)
    for adjusted in [True, False]:
        record, points = fit_model(table, adjusted = adjusted)
        records.append(record)
        if adjusted:
            np.testing.assert_allclose(
                record['standardizedCoefficient'], -0.24064684758478283
            )
            assert round(record['pValueUnadjusted'], 4) == 0.0118
            points.to_csv(
                args.output / 'Figure_3D_source_data.csv', index = False
            )
    result = pd.DataFrame(records)
    result.to_csv(
        args.output / 'matched_utr3_models.tsv', sep = '\t', index = False
    )
    print(result.to_string(index = False))


if __name__ == '__main__':
    main()
