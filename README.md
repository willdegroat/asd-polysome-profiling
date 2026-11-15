# ASD polysome profiling

Analysis source code for paired input, monosome-associated, and
polysome-associated RNA sequencing from NGN2-induced neurons. The study contains
one individual with ASD and an unaffected sibling, with three independently
processed neuronal cultures per individual. Each culture contributes matched
Input, Monosome, and Polysome libraries. Culture replicates are not independent
donors, and these data do not support population-wide autism effects.

## Published results

The frozen numerical results are included as compressed files in `published/`.
Restore and verify them without refitting the models:

```console
python scripts/restore_results.py
python scripts/verify_manuscript.py
```

See [manuscript results](manuscript/README.md) for the individual figure panels,
Table 1, source data, and supplementary analyses. Render panels separately;
the repository does not assemble them into full figures.

## Workflow

Run the notebooks in numeric order. To execute the complete workflow and save
the executed notebooks in place:

```bash
python scripts/run_notebooks.py
```

| Notebook | Analysis |
| --- | --- |
| `00_data.ipynb` | Metadata and nf-core/rnaseq provenance |
| `01_counts.ipynb` | Salmon matrix validation |
| `02_qc.ipynb` | Library and sample-level quality control |
| `03_expression.ipynb` | Input RNA differential expression |
| `04_fractions.ipynb` | Monosome and polysome differential abundance |
| `05_redistribution.ipynb` | Paired polysome–monosome modeling with Monopoly |
| `06_sensitivity.ipynb` | Culture, normalization, and label sensitivity |
| `07_enrichment.ipynb` | Ranked enrichment and overrepresentation analysis |
| `08_integration.ipynb` | Cross-layer gene and pathway integration |
| `09_transcripts.ipynb` | Differential transcript abundance |
| `10_usage.ipynb` | Differential transcript usage with satuRn and stageR |
| `11_architecture.ipynb` | Transcript and 5′ UTR architecture |
| `12_convergence.ipynb` | Fu et al. ASD-gene convergence |

Gene-level models use Salmon length-scaled counts and PyDESeq2. Monopoly fits
the participant-paired condition-by-fraction model. Transcript-abundance models
reconstruct tximport-compatible `lengthScaledTPM` counts, while transcript-usage
models reconstruct `dtuScaledTPM` counts before satuRn and stageR inference.

## Environment

Python 3.12 was used for the analysis.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

[Monopoly 0.1.0](https://github.com/willdegroat/monopoly) must be installed in
the same environment. For a sibling clone:

```bash
python -m pip install -e ../monopoly
```

Notebook 10 additionally requires R 4.6, Bioconductor 3.23, satuRn 1.20.0,
stageR 1.34.0, and BiocParallel 1.46.0. Exact R session information is written
with the DTU results.
`Rscript` must be available on `PATH`.

## Data and results

The required input layout is documented in [`data/README.md`](data/README.md).
Large quantification and reference files are intentionally excluded from Git.
Pinned gene-set collections and deidentified sample metadata are versioned.

All derived tables and diagnostic files are written to `results/`. Frozen copies
in `published/` retain all analyses, including null findings and sensitivity
checks. Large quantification matrices and reference sequences are not included;
the required input files and their SHA-256 fingerprints are documented separately.

Monopoly uses a shared-dispersion working model. The benchmark demonstrates that
unequal dispersions can inflate false-discovery rates. RNA redistribution does
not directly measure protein synthesis or ribosome occupancy. See
[interpretation and reproducibility](manuscript/README.md).

## Authors

William DeGroat and Anat Kreimer.

## License

This project is released under the MIT License.

## References

- DESeq2: <https://doi.org/10.1186/s13059-014-0550-8>
- PyDESeq2: <https://doi.org/10.1093/bioinformatics/btad547>
- GSEA: <https://doi.org/10.1073/pnas.0506580102>
- SynGO: <https://doi.org/10.1016/j.neuron.2019.05.002>
- Fu et al. ASD risk genes: <https://doi.org/10.1038/s41588-022-01104-0>
- tximport: <https://doi.org/10.12688/f1000research.7563.1>
- satuRn: <https://doi.org/10.12688/f1000research.51749.2>
- stageR: <https://doi.org/10.1186/s13059-017-1277-0>
