"""Verify the numerical findings reported in the manuscript figures and table."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, norm

ROOT = Path(__file__).resolve().parents[1]


def read_result(name: str) -> pd.DataFrame:
    """Read a restored analysis table."""
    return pd.read_csv(ROOT / 'results' / name, sep = '\t')


def verify() -> None:
    """Check frozen files, discovery counts, selections, and testing families."""
    manifest = json.loads((ROOT / 'published/manifest.json').read_text())
    for name, record in manifest['inputs'].items():
        path = ROOT / name
        if path.exists():
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            assert digest == record['sha256'], path
    for record in manifest['files']:
        path = ROOT / record['path']
        if not path.exists():
            raise FileNotFoundError(
                'Run python scripts/restore_results.py first.'
            )
        assert (
            hashlib.sha256(path.read_bytes()).hexdigest() == record['sha256']
        ), path

    metadata = pd.read_csv(ROOT / 'data/metadata.tsv', sep = '\t')
    assert len(metadata) == 6
    assert metadata.groupby('condition').size().eq(3).all()
    input_data = read_result('inputDifferentialExpression.tsv')
    pm = read_result('polysomeMonosome.tsv')
    for table, expected in [(input_data, (135, 170)), (pm, (12, 10))]:
        significant = table[table.adjustedPValue.lt(0.05)]
        counts = (
            significant.log2FoldChange.gt(0).sum(),
            significant.log2FoldChange.lt(0).sum(),
        )
        assert counts == expected, counts
        finite = table.pValue.notna() & table.statistic.notna()
        np.testing.assert_allclose(
            table.loc[finite, 'pValue'],
            2 * norm.sf(abs(table.loc[finite, 'statistic'])),
            atol = 1e-14,
        )

    ora = read_result('overrepresentationAnalysis.tsv')
    go = ['goBiologicalProcess', 'goCellularComponent', 'goMolecularFunction']
    selected = ora[
        ora.analysis.eq('input')
        & ora.collection.isin(go)
        & ora.adjustedPValue.lt(0.05)
    ].sort_values(
        ['adjustedPValue', 'pValue', 'collection', 'term', 'direction']
    )
    assert selected.head(6).direction.eq('negative').all()
    assert len(selected.head(6)) == 6

    annotations = json.loads(
        (ROOT / 'manuscript/assets/gene_biotypes.json').read_text()
    )['records']
    coding = pm[
        pm.adjustedPValue.lt(0.05)
        & pm.geneId.map(
            lambda gene: annotations.get(gene, {}).get('geneType')
        ).eq('protein_coding')
    ].copy()
    coding['absoluteEffect'] = coding.log2FoldChange.abs()
    coding = coding.sort_values(
        ['absoluteEffect', 'geneId'], ascending = [False, True]
    )
    assert coding.head(10).geneName.tolist() == [
        'H2BC17',
        'DHRSX',
        'DNAI3',
        'TOP2A',
        'NPIPB3',
        'PHOX2A',
        'NEUROG1',
        'SLC31A1',
        'GCN1',
        'ATP1B2',
    ]
    gene_effects = pm.set_index('geneName')
    assert round(2 ** gene_effects.loc['PHOX2A', 'log2FoldChange'], 1) == 2.9
    assert round(2 ** -gene_effects.loc['TOP2A', 'log2FoldChange'], 1) == 4.8

    focused = read_result('fuAsdRankedEnrichment.tsv')
    assert len(focused) == 8
    np.testing.assert_allclose(
        focused.globalAdjustedPValue,
        false_discovery_control(focused.pValue),
    )
    omission = read_result('fuAsdLeaveOneOutEnrichment.tsv')
    for _, table in omission.groupby('tier'):
        assert len(table) == 6
        np.testing.assert_allclose(
            table.tierAdjustedPValue,
            false_discovery_control(table.pValue),
        )
        missing = table.loc[
            table.tierAdjustedPValue.ge(0.05), 'omittedParticipantId'
        ]
        assert missing.astype(str).tolist() == ['1851']
    enrichment = read_result('rankedEnrichment.tsv')
    significant = enrichment[
        enrichment.analysis.eq('redistribution')
        & enrichment.collection.isin(['synGO', 'fuAsd'])
        & enrichment.adjustedPValue.lt(0.05)
    ]
    assert len(significant) == 9

    evidence = read_result('transcriptGeneEvidence.tsv')
    evidence = evidence[
        evidence.analysis.eq('redistribution') & evidence.variant.eq('stable')
    ]
    assert len(evidence) == 11971
    assert (evidence.geneFdr & ~evidence.parentGeneFdr).sum() == 48
    assert (evidence.geneFdr & evidence.parentGeneFdr).sum() == 5
    assert (~evidence.geneFdr & evidence.parentGeneFdr).sum() == 9
    np.testing.assert_allclose(
        evidence.geneAdjustedPValue,
        false_discovery_control(evidence.genePValue),
    )
    features = read_result('transcriptFeatureAssociation.tsv')
    assert len(features) == 28
    np.testing.assert_allclose(
        features.adjustedPValue,
        false_discovery_control(features.pValue),
    )
    pm_features = features[features.analysis.eq('redistribution')]
    assert pm_features.adjustedPValue.lt(0.05).sum() == 5
    assert pm_features.transcripts.eq(11679).all()
    assert pm_features.genes.eq(3843).all()

    for scenario in ['count_model_grid', 'equal_n20']:
        folder = ROOT / 'manuscript/benchmarks' / scenario
        rows = pd.read_csv(folder / 'replicateResults.tsv', sep = '\t')
        if scenario == 'count_model_grid':
            rows = rows[rows.scenarioId.eq('unequalEffect1N20')]
        assert len(rows) == 800 and rows.methodError.isna().all()
        assert rows.groupby('method').size().eq(200).all()
    print('All frozen files and manuscript numerical checks passed.')


if __name__ == '__main__':
    verify()
