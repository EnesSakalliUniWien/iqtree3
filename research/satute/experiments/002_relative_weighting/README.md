# 002 - Relative eigenvalue weighting

Runs paired comparisons between the published dominant statistic and relative
eigenvalue weighting under JC, GTR, and LG settings.

- Driver: `run.py`
- Shard combiner: `combine.py`
- Renderers: `render.py`, `render_power_summary.py`
- Config: `../../config/experiments/002_relative_weighting.yaml`

The production driver compares the native IQ-TREE `dominant` and
`eigenvalue_weighted` rows directly. It does not recompute either statistic in
Python for every replicate; the independent Python implementation is retained
only in `../../tests/native/` for validation. The production path neither loads
nor imports NumPy. The replicate seed is used for both simulation and IQ-TREE
analysis, so topology-search reruns are deterministic. Shard detail rows are
flushed after every replicate task so `--resume` can recover completed work
after interruption.
