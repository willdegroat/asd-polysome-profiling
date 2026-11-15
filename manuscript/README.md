# Manuscript results

The main results comprise Figures 1â€“3 and Table 1. All comparisons refer to one
ASD individual and an unaffected sibling, with three cultures per background.

## Reproduce the publication outputs

Install the analysis environment and Monopoly, then restore the frozen tables:

```console
python -m pip install -r requirements.txt
python scripts/restore_results.py
python scripts/verify_manuscript.py
python scripts/render_figures.py
python scripts/check_figures.py
python scripts/make_table_1.py
```

Each figure panel is saved separately under `figures/generated/`. Select panels
with `--panels 2A 2B`, or choose an output directory with `--output`. Arial is
required to reproduce the approved typography. The renderer stops if Arial is
unavailable rather than silently substituting another font.

The publication source data are under `manuscript/source_data/`. The processed
RNA-sequencing quantifications are restored from the archive. Large GENCODE
references are needed for sequence-feature reconstruction, but not for inspecting
the archived numerical results.
Input file fingerprints are included in `published/manifest.json`.

## Main analyses

| Output | Analysis |
| --- | --- |
| Figure 1 | Study overview and fraction-specific and combined PCA |
| Figure 2Aâ€“B | Input RNA abundance and GO overrepresentation |
| Figure 2Câ€“F | Gene-level P/M effects and PHOX2A/TOP2A expression profiles |
| Figure 2Gâ€“I | Fu ASD gene-set and SynGO enrichment, with culture omissions |
| Figure 3Aâ€“B | Transcript-derived gene evidence and the DLK1 example |
| Figure 3C | Joint within-gene sequence-architecture model |
| Figure 3D | Exact-CDS/5â€²-UTR-matched 3â€² UTR association |
| Table 1 | Equal- and unequal-dispersion benchmarks at 20 per condition |

Figure 2B uses the six lowest adjusted P values across GO collections and both
directions, without manually adding neural terms. Figure 2C selects the ten
significant protein-coding genes with the largest absolute P/M effects.
Expression profiles average logâ‚‚(TPM + 1), not raw TPM. Their P/M annotations
come from the fitted count model, not ratios of the displayed means.

Figure 2G uses BH correction across eight focused Ensembl-ID enrichment tests.
Figure 2H corrects across six omissions within each gene set. Figure 2I uses
the general gene-symbol workflow and collection-specific GSEA FDR q values.
These are different analysis populations and testing families.

Figure 3A counts genes with both tests available. Its 14 gene-level discoveries
therefore differ from the 22 discoveries in the complete gene-level analysis.
Transcript-derived gene evidence is not a formal test of isoform switching.
Figure 3C retains the original correction across 28 feature-by-layer tests,
although only the seven P/M associations are displayed.

Figure 3D compares 154 transcripts from 69 genes in 70 groups with identical
coding and 5â€² UTR sequences. The model centers variables within matching groups,
adjusts for transcript abundance and mean within-gene transcript proportion, and
uses parent-gene-clustered standard errors. Its outcome is P/M logâ‚‚ fold change;
Figure 3C instead models standardized Wald statistics. The displayed follow-up
P value is unadjusted. This association is not significant without the abundance
and proportion covariates. The follow-up scripts and all fitted alternatives are
retained for transparency.

## Supplementary results and limitations

The archive includes all fraction abundance analyses, gene and transcript
redistribution, GO/SynGO/Fu enrichment, transcript usage (including null results),
normalization and culture sensitivity, label perturbations, and sequence-feature
analyses. The benchmark archive includes the complete original count-model grid,
FracFixR native-estimand comparisons, and the additional equal-dispersion run.

ASD gene-set enrichment loses significance when culture ASD1 is omitted.
Controlled culture-omission refits also show sensitivity of the architecture
associations, including a change in the broad 3â€² UTR association after omitting
U3. These checks do not justify excluding a culture or claiming a causal sequence
effect. Negative Input tests do not establish unchanged RNA abundance.

Monopoly reports relative association with the sequenced pools, not absolute RNA
partitioning or protein production. The Monosome pool includes 40S, 60S, and 80S
regions and is not assumed to be translationally inactive. Inference uses a
shared-dispersion working model; unequal-dispersion benchmarks show FDR inflation.
Table 1 reports mean per-dataset false-discovery proportion. Pooled FDR and Monte
Carlo standard errors remain in the complete supporting benchmark results.

Fraction-specific PCA uses R DESeq2 1.52.0 blind VST followed by centered PCA on
all nonzero genes. The combined PCA uses logâ‚‚(TPM + 1) and excludes zero-variance
genes. No differential-expression selection is used for PCA.

To recompute each final VST from the raw matrices, use
`Rscript scripts/prepare_pca.R Input` and repeat for `Monosome` and `Polysome`.
These checks write to `manuscript/work/pca/` and leave the archived figure inputs
unchanged. The numbered QC notebook uses the Python implementation; the final
publication panels use the explicitly archived R DESeq2 transformation instead.
