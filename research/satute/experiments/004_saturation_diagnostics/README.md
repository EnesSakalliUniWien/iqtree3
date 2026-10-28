# 004 - Saturation diagnostics

Separates raw coherence, internal uncertainty, standardized detection,
category persistence, and the monotonic model-based saturation index for the
two production formulas (`dominant` and `eigenvalue_weighted`).

Pooled rows use the integrated mixture score and its sample variance.
Hard-category rows are diagnostics with no nominal p-values. The historical
`mean_pooled_contribution` column retains its schema name and records
`(valid sites / all sites) × category satC`; summing that diagnostic column
does not reconstruct the integrated pooled score.

- Driver: `run.py`
- Null summary renderer: `summarize_null.py`
- Manifest and run configuration: `experiment.json`
- Common workflow adapter: `workflow.py`
