# 002 - Relative eigenvalue weighting

Runs paired comparisons between the published dominant statistic and relative
eigenvalue weighting under nucleotide and protein settings.

- Driver: `run.py`
- Shard combiner: `combine.py`
- Renderers: `render.py`, `render_power_summary.py`,
  `render_formula_comparison.R`
- Config: `../../config/experiments/002_relative_weighting.yaml`

The production detail schema is versioned and keeps the three native decision
contracts side by side: `decision_unadjusted`, `decision_taxon_bonf`, and
`decision_fdr`. The normalized summary adds a `decision_rule` dimension, so a
Bonferroni or BY-FDR plot never requires a second simulation. Every shard also
writes `head_to_head_branch_audit.tsv`, `head_to_head_fdr_audit.tsv`, and
`head_to_head_timing.tsv`; the first is required for empirical tree-wide FDR
audits and the last identifies simulation, topology-search, and I/O
bottlenecks. The audit files alone do not estimate empirical FDR: a publication
claim about FDR control still requires a separate truth-labelled exact-null and
mixed-null calibration design, with held-out FDP/FDR/TPR summaries. Until that
design is run, `by_fdr` is reported as a native BY decision rule, not as a
validated operating-characteristic claim.

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

## Independent 2D and 3D formula plots

`render_formula_comparison.R` produces three independent 2D figures and three
independent 3D line plots:

- published dominant coefficient;
- eigenvalue-weighted coefficient;
- paired difference (`eigenvalue_weighted - dominant`).

The 2D figures are written as vector PDF/SVG and 450-DPI PNG. Each 3D figure is
written both as a multi-panel vector PDF containing every tree/scenario view
and as a self-contained interactive Plotly HTML widget. Target branch length,
alignment length, and informative fraction are the three axes. The figures draw
only measured alignment-length trajectories and do not interpolate a surface.
The difference figures pair formula decisions on the same simulated alignment
and report a 95% paired bootstrap interval plus the exact McNemar p-value in a
companion TSV.

```bash
Rscript research/satute/experiments/002_relative_weighting/render_formula_comparison.R \
  --summary /path/to/merged/head_to_head_summary.tsv \
  --detail /path/to/merged/head_to_head_detail.tsv \
  --outdir /path/to/figures \
  --simulation-model GTR_PF06346 \
  --evaluation-model GTR_PF06346 \
  --expected-reps 250
```

Pass `--decision-rule unadjusted`, `--decision-rule taxon_bonferroni`, or
`--decision-rule by_fdr` to render the corresponding correction contract.
`by_fdr` is the tree-wide Benjamini--Yekutieli result; it is not a relabeling
of the paper's taxon-pair Bonferroni correction.

The separate calibration gate and its 5,000 exact-null / 2,000 mixed-null
design are specified in
[`docs/experiments/fdr-calibration-design.md`](../../docs/experiments/fdr-calibration-design.md).

The publication curve design uses fixed aliases rather than fitting a rate
mixture independently in every replicate. `GTR_PF06346_G4` resolves to the
frozen PF06346 exchangeabilities and frequencies with `+G4{0.5}`;
`GTR_PF06346_I_G4` adds the fixed invariant proportion `+I{0.1}`. The latter
must use IQ-TREE AliSim. The driver rejects it under Seq-Gen instead of silently
dropping the invariant component.

Truth-labelled multiple-testing calibration is run separately with
`fdr_calibration_run.py`. It simulates exact independence across the target
edge under homogeneous fixed JC/GTR models, retains normalized-split truth for
every pooled branch, and reports per-replicate FDP and power for unadjusted,
taxon-Bonferroni, and BY rules. The 5,000-replicate exact-null and
2,000-replicate mixed-null designs are submitted with
`tools/cluster/lisc/submit_fdr_calibration_suite.sh`; their outputs must not be
merged into the publication power-curve summaries.

The publication-summary renderer also requires schema-v2 summaries and accepts
the same `--decision-rule` flag. It refuses incomplete replicate cells unless
`--allow-incomplete` is explicitly supplied, and it derives the observed site,
branch and tree-case grid instead of silently filling missing cells from the
legacy Fig. 2 defaults.

For regenerated rate-heterogeneous curves, pass separate summaries and exact
model strings with `--gtr-summary`, `--lg-summary`,
`--gtr-simulation-model`, `--gtr-evaluation-model`, `--lg-simulation-model`,
and `--lg-evaluation-model`; the legacy `--rate-summary` option remains an
alias for the GTR summary.

Required R packages are `ggplot2`, `plot3D`, `plotly`, `ragg`, `svglite`,
`scales`, `viridisLite`, and `htmlwidgets`.
