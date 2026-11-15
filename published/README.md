# Frozen results

Each gzip file preserves the original numerical output byte for byte.
`manifest.json` records uncompressed and compressed SHA-256 checksums, file sizes,
and the upstream input fingerprints. Restore the files with:

```console
python scripts/restore_results.py
python scripts/verify_manuscript.py
```

Restoration refuses to overwrite results that differ from the archive unless
`--overwrite` is explicitly supplied. Rerunning models can change diagnostic
metadata, runtime, or last-bit numerical values; the frozen archive remains the
record used by the manuscript.

The archive includes all primary analysis outputs, supplementary/null findings,
architecture follow-up models and controlled refits, individual-panel source data,
and all three benchmark result sets. It also includes all processed quantification matrices, MultiQC reports, and
the Fu supplementary workbook. It excludes raw sequence reads, large reference
files, local Python environments, and temporary renderings.

The original combined-PCA export had a duplicate `sampleId` header. The untouched
export is retained, alongside `manuscript/source_data/1E/PCA_scores.csv`, which
removes the row counter and gives the culture identifier one unambiguous header.
No coordinates or labels were changed.
