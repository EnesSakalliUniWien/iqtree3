# 002 - Relative eigenvalue weighting

Runs paired comparisons between the published dominant statistic and relative
eigenvalue weighting under JC, GTR, and LG settings.

- Driver: `run.py`
- Shard combiner: `combine.py`
- Renderers: `render.py`, `render_power_summary.py`
- Config: `../../config/experiments/002_relative_weighting.yaml`

The driver uses the replicate seed for both simulation and IQ-TREE analysis,
so topology-search reruns are deterministic. Its independent Python reference
compresses identical alignment patterns with exact frequency weighting and
reuses the inferred-tree calculation for unadjusted and Bonferroni decisions.
Shard detail rows are flushed after every replicate task so `--resume` can
recover completed work after interruption.
