"""Create Table 1 from the complete frozen benchmark summaries."""

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
METHODS = ['monopoly', 'edgeR', 'ratioWelch']
NAMES = {
    'monopoly': 'Monopoly',
    'edgeR': 'edgeR QL',
    'ratioWelch': 'Ratio Welch test',
}


def main() -> None:
    """Export the displayed values and their unrounded source rows."""
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument('--output', type = Path, default = ROOT / 'manuscript/work')
    args = parser.parse_args()
    args.output.mkdir(parents = True, exist_ok = True)
    rows = []
    for label, folder, scenario in [
        ('Equal', 'equal_n20', 'equalEffect1N20'),
        ('Unequal', 'count_model_grid', 'unequalEffect1N20'),
    ]:
        path = ROOT / 'manuscript/benchmarks' / folder / 'metricSummary.tsv'
        table = pd.read_csv(path, sep = '\t')
        table = (
            table[table.scenarioId.eq(scenario)]
            .set_index('method')
            .loc[METHODS]
        )
        assert (
            table.replicates.eq(200).all()
            and table.failedReplicates.eq(0).all()
        )
        table = table.reset_index()
        table.insert(0, 'Dispersion', label)
        rows.append(table)
    source = pd.concat(rows, ignore_index = True)
    source.to_csv(
        args.output / 'Table_1_source_data.tsv', sep = '\t', index = False
    )
    displayed = []
    for row in source.itertuples():
        runtime = (
            f'{row.runtimeSecondsMean:.3f}'
            if row.runtimeSecondsMean < 0.1
            else f'{row.runtimeSecondsMean:.1f}'
        )
        displayed.append(
            [
                row.Dispersion,
                NAMES[row.method],
                f'{100 * row.powerMean:.1f}',
                f'{100 * row.falseDiscoveryRateMean:.1f}',
                f'{row.effectRmseMean:.3f}',
                runtime,
            ]
        )
    headers = [
        'Dispersion',
        'Method',
        'Power (%)',
        'FDR (%)',
        'Effect RMSE',
        'Runtime (s)',
    ]
    pd.DataFrame(displayed, columns = headers).to_csv(
        args.output / 'Table_1.tsv',
        sep = '\t',
        index = False,
    )
    text = '| ' + ' | '.join(headers) + ' |\n|---|---|---:|---:|---:|---:|\n'
    text += '\n'.join('| ' + ' | '.join(row) + ' |' for row in displayed) + '\n'
    (args.output / 'Table_1.md').write_text(text, encoding = 'utf-8')
    print(text)


if __name__ == '__main__':
    main()
