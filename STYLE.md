# Analysis style

Keep the analysis explicit and readable. Use single-quoted Python strings,
88-character lines, spaces around keyword and default assignments, and concise
comments explaining scientific choices. Use snake_case internally and camelCase
in exported scientific tables. Preserve complete testing families, identifiers,
null results, and sensitivity analyses. Do not silently select favorable results.

Use the existing notebooks for the full analysis, plain scripts for publication
outputs, and established numerical libraries for model fitting. Figure scripts
render individual panels; they must not assemble a full manuscript figure.
