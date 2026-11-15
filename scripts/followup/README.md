# Architecture follow-up

The frozen outputs are restored into `manuscript/followup/`. These analyses
include all investigated covariates and matched populations, rather than only
the displayed association. Reproduce Figure 3D directly with:

```console
python scripts/matched_utr3.py
```

To rerun the broader investigation, first provide the raw quantifications and
run-matched transcript FASTA described in `data/README.md`. Run from the repository
root, in this order:

```console
python scripts/followup/investigate_utr3.py
python scripts/followup/investigate_utr3_stage2.py
python scripts/followup/investigate_utr3_stage3.py
python scripts/followup/reconstruct_count_model.py
Rscript scripts/followup/refit_fixed_dispersion_NB.R
python scripts/followup/analyze_all_NB_omissions.py
Rscript scripts/followup/model_inference.R
```

The R refits require DESeq2 1.52.0. They retain the original size factors and
dispersions and are controlled sensitivity analyses, not exact reruns of the
PyDESeq2 primary estimates. Completed R refits are reused; move those seven
derived CSV files aside before intentionally recomputing them. The scripts write
derived outputs to `manuscript/followup/`, leaving `published/` unchanged.

The source sequence cache is validated against the reference file identity and
transcript population before reuse. The current working directory must be the
repository root for the R scripts.
