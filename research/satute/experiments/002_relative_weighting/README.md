# 002 - Relative eigenvalue weighting

Runs paired comparisons between the published dominant statistic and relative
eigenvalue weighting under nucleotide and protein settings.

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

The fixed protein extension provides `LG`, `WAG`, `JTT`, and `Q.PFAM`
aliases, plus `LG_G4`, `WAG_G4`, `JTT_G4`, and `Q.PFAM_G4`. The `*_G4`
aliases use the explicitly fixed shape `+G4{0.5}` so matched simulations do
not estimate rate heterogeneity separately in every replicate. Use IQ-TREE
AliSim for the complete extension because Seq-Gen does not implement Q.PFAM.
A `GTR20_REF` alias is added only after its 189 exchangeabilities, 20
frequencies, and gamma shape have been estimated once from the independent
reference alignment and recorded with checksums.
