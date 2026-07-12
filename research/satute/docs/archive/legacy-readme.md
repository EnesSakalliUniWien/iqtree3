# SatuTe Wiki

This directory documents the SatuTe paper, the IQ-TREE integration work in
this repository and the proposed eigenvalue-weighted extension. It separates
source material, derivations, reproducible development checks and manuscript
drafts. The reduced simulation files are used to verify software behavior; they
are not evidence for the empirical claims made by the SatuTe paper.

## Source Material

The paper is:

Manuel C, Sakalli E, Schmidt HA, Viñas C, von Haeseler A, Elgert C. 2025.
When the Past Fades: Detecting Phylogenetic Signal with SatuTe. Molecular
Biology and Evolution 42(5):msaf090. https://doi.org/10.1093/molbev/msaf090

Local copy:

```text
doc/satute-wiki/sources/msaf090.pdf
doc/satute-wiki/sources/SatuTe_manuscript_2024_08/
```

Online sources checked:

- Oxford Academic article: https://academic.oup.com/mbe/article/42/5/msaf090/8151464
- SatuTe implementation: https://github.com/Elli-ellgard/SatuTe
- SatuTe example analyses: https://github.com/Elli-ellgard/SatuTe-example-analyses
- EvoNAPS web interface: http://evonaps.cibiv.univie.ac.at/

## Versioned Artifact Policy

This directory versions source notes, scripts, the addendum paper source/PDF,
small checked development-run summaries, original reference figures used by the
addendum, and recovered public source material. It does not version raw
simulation run directories, generated full-run figures, local LaTeX build
products, or vendored Seq-Gen source/binaries. Those paths are ignored by
`doc/satute-wiki/.gitignore`; regenerate them from the scripts and pass
`--seqgen /path/to/seq-gen` for paper-faithful Seq-Gen reruns.

## Reading Order

1. [01-paper-map.md](01-paper-map.md) summarizes the paper, code, data availability,
   and simulation design.
2. [02-formula-derivation.md](02-formula-derivation.md) derives the
   posterior/eigenvector statistic from the branch-side likelihood vectors.
3. [03-reconstruct-paper-simulations.md](03-reconstruct-paper-simulations.md)
   documents the full Fig. 2 simulation target and the reduced development
   script used for command-line verification.
4. [04-new-eigenvalue-statistic.md](04-new-eigenvalue-statistic.md) defines
   the eigenvalue-weighted statistic as an extension of the published
   dominant-mode statistic.
5. [05-saturation-diagnostics.md](05-saturation-diagnostics.md) audits the
   native implementation choices and separates coherence, uncertainty,
   detection, category persistence and a monotonic model-based saturation
   index.
6. [06-rotation-relationship.md](06-rotation-relationship.md) relates the
   monotonic scale to the spectral angles and branch construction in *When
   Reversible Mutation Becomes Rotation*, and maps the dominant, relative
   weighted, soft-mixture and model-based quantities to their distinct roles.
7. [saturation_diagnostics.py](saturation_diagnostics.py) runs the focused
   fixed-model validation and generates its diagnostic figures and assumption
   report. Native `.sat.stat` output includes the formula-independent
   `InformationFraction` and monotonic `SaturationIndex` columns, with pooled
   `satInfo` and `satIndex` annotations in the NEXUS tree.
8. [verify_saturation_null.py](verify_saturation_null.py) enumerates the exact
   saturated pattern distribution and verifies null centering, scale and type-I
   error across alignment lengths.
9. [verify_weighted_likelihood.py](verify_weighted_likelihood.py) verifies the
   optimal persistence weights, the positive-rate soft-mixture likelihood
   identity, and the corrected invariant-category null by exact enumeration.
10. [weighted_invariant_simulations.py](weighted_invariant_simulations.py)
   runs the dedicated paired `+I+G4`/`+I+R4` exact-pattern grid, including
   held-out null calibration and paired bootstrap power intervals.
11. [render_weighted_invariant_shift.py](render_weighted_invariant_shift.py)
   plots the direct paired shifts in the null statistic, false-positive rate,
   calibrated decision margin and power.
12. [render_weighted_invariant_full.py](render_weighted_invariant_full.py)
   renders the dense 27-step branch analysis as publication-quality calibrated
   power curves, paired-difference heatmaps and null-calibration panels in
   vector PDF/SVG and 3x PNG formats.
13. [head_to_head_satute_simulations.py](head_to_head_satute_simulations.py)
   is the paired simulation driver: one simulated alignment is evaluated by
   the published dominant and eigenvalue-weighted tests under the same scenario
   and threshold.
14. [short-addendum-paper.tex](short-addendum-paper.tex) is a manuscript draft
   for the eigenvalue-weighted extension. It uses the SatuTe notation and does
   not present the reduced development checks as empirical evidence. The
   compiled PDF is [short-addendum-paper.pdf](short-addendum-paper.pdf).
15. [five-page-explainer.md](five-page-explainer.md) is the longer working
   explainer kept for background.

## Run The Paired Head-To-Head Smoke Test

Use a Python runtime with NumPy installed:

```bash
${PYTHON_BIN:-python3} \
  doc/satute-wiki/head_to_head_satute_simulations.py \
  --iqtree ./build/iqtree3 \
  --outdir /tmp/iqtree-satute-head-to-head \
  --reps 1 \
  --site-lengths 100 \
  --branch-lengths 0.1 \
  --tree-cases five_external,sixteen_internal \
  --simulation-models JC \
  --evaluation-models JC \
  --scenario-set fig2 \
  --simulator alisim
```

The output has one row per replicate, analysis scenario and formula. The same
alignment and target branch are evaluated by the published dominant and
eigenvalue-weighted statistics.

The verified development-check outputs are:

```text
doc/satute-wiki/results/head_to_head_development_detail.tsv
doc/satute-wiki/results/head_to_head_development_summary.tsv
```

## Biological Rerun Results

The curated biological rerun is:

```text
doc/satute-wiki/results/biological/enhanced_sliding_window_rerun_20260619_150104/
```

This run restores the 16S rRNA Tree-of-Life sliding-window analysis for the
branch leading to Eukaryota and the branch leading to yeast. It contains the
published SatuTe rerun, the native dominant-mode calculation and the
eigenvalue-weighted calculation, with the comparison figures:

```text
comparison_to_manuscript/original_vs_enhanced_sliding_window_curves.pdf
comparison_to_manuscript/original_vs_enhanced_saturated_windows.pdf
comparison_to_manuscript/original_vs_enhanced_summary.tsv
```

The rerun uses the manuscript GTR+F+G4 nucleotide model and four rate
categories. The protein LG+G4 Tree-of-Life dataset was not recomputed with the
enhanced statistic in this run.

The comparison figures can be regenerated from the curated tables with:

```bash
python3 doc/satute-wiki/plot_biological_comparison.py
```

The independent biological comparison is:

```text
doc/satute-wiki/results/biological/independent_enhanced_sliding_window_20260619_151135/
```

This run includes both manuscript Tree-of-Life examples. The protein analysis
uses `protein_based_2D_tree` with `trees/protein_based_2D_ToL.treefile`; the
16S rRNA analysis uses `rRNA_based_3D_tree` with
`trees/rRNA_based_3D_ToL.treefile`. In the protein 2D tree, the branch leading
to Eukaryota has 2,561 windows, with 138 saturated windows in the published
SatuTe rerun and 5 in the eigenvalue-weighted calculation. The yeast branch has
2,561 windows, with counts 54 and 8. In the 16S rRNA 3D tree, the branch
leading to Eukaryota has 1,912 windows, with counts 822 and 699. The yeast
branch has 1,912 windows, with counts 70 and 29.

The independent comparison figures can be regenerated with:

```bash
python3 doc/satute-wiki/plot_independent_biological_comparison.py
```

## Paper-Faithful Rerun Target

The original manuscript archive specifies Seq-Gen v1.3.4, JC simulation,
1,000 alignments per point, sequence lengths `100,1000,10000`, the 16-value
branch grid, and four paired analysis scenarios. The five-taxon topology is a
single taxon `B` attached to the root of a balanced four-taxon `T_A` with all
subtree branches set to `0.2`. The 16-taxon topology uses balanced
`A1..A8 | B1..B8` subtrees with internal and external subtree branch lengths
sampled from EvoNAPS distributions.

Exact rerun command shape:

```bash
${PYTHON_BIN:-python3} \
  doc/satute-wiki/head_to_head_satute_simulations.py \
  --iqtree ./build/iqtree3 \
  --seqgen /path/to/seq-gen \
  --simulator seq-gen \
  --evonaps-branch-lengths /path/to/evonaps_branch_lengths.tsv \
  --outdir /tmp/iqtree-satute-head-to-head-full \
  --reps 1000 \
  --paper-grid \
  --tree-cases five_external,sixteen_internal \
  --simulation-models JC \
  --evaluation-models JC \
  --scenario-set fig2
```

This command is JC-only because Fig. 2 in the paper simulated and evaluated
under JC. Do not include a GTR model in the paper-faithful Fig. 2 rerun.

The paper's GTR experiment is a separate Supplementary Fig. model
misspecification study. Its simulation model is available as
`GTR_PF06346`, corresponding to the EvoNAPS PF06346 best-fit model:

```text
rates AC=0.6676, AG=3.7807, AT=4.2833, CG=0.5354, CT=0.8718, GT=1.0
frequencies A=0.125, C=0.436, G=0.191, T=0.245
```

That supplementary experiment should simulate with `GTR_PF06346` and evaluate
separately under the correctly specified GTR model and the misspecified JC,
K2P and F81 models. It should not be mixed into the Fig. 2 reproduction.

Supplementary command shape:

```bash
${PYTHON_BIN:-python3} \
  doc/satute-wiki/head_to_head_satute_simulations.py \
  --iqtree ./build/iqtree3 \
  --seqgen /path/to/seq-gen \
  --simulator seq-gen \
  --evonaps-branch-lengths /path/to/evonaps_branch_lengths.tsv \
  --outdir /tmp/iqtree-satute-misspecification-full \
  --reps 1000 \
  --paper-grid \
  --tree-cases five_external,sixteen_internal \
  --simulation-models GTR_PF06346 \
  --evaluation-models GTR_PF06346,JC,K2P,F81 \
  --scenario-set misspecification
```

Skewed GTR analyses are extension experiments, not reproductions of the
published simulation figures. Run them as explicit model pairs so that each
simulation model is evaluated under the intended model rather than as an
accidental cross-product:

```bash
${PYTHON_BIN:-python3} \
  doc/satute-wiki/head_to_head_satute_simulations.py \
  --iqtree ./build/iqtree3 \
  --seqgen /path/to/seq-gen \
  --simulator seq-gen \
  --evonaps-branch-lengths /path/to/evonaps_branch_lengths.tsv \
  --outdir /tmp/iqtree-satute-skewed-gtr-extension \
  --reps 1000 \
  --paper-grid \
  --tree-cases five_external,sixteen_internal \
  --model-pairs GTR_SKEW_FREQ:GTR_SKEW_FREQ,GTR_SKEW_RATES:GTR_SKEW_RATES,GTR_SKEW_BOTH:GTR_SKEW_BOTH \
  --scenario-set all
```

The EvoNAPS table must contain a branch type column named `type`, `kind`, or
`branch_type` with values starting with `internal` or `external`, and a length
column named `length`, `branch_length`, `blen`, or `bl`. The official EvoNAPS
download uses `BRANCH_TYPE` values `i/e` and `BL`, which the driver accepts.

Recovered EvoNAPS public branch source:

```text
doc/satute-wiki/sources/evonaps/evonaps_16taxon_branch_lengths.tsv
```

This TSV was downloaded from the EvoNAPS tree search result page for DNA trees
with 16 leaves and maximum internal/external branch lengths no greater than
`0.2`. After applying the manuscript interval filter (`0.04 <= length <= 0.2`),
the usable pool contains 92 internal and 233 external branch lengths. Attempts
to download the unfiltered all-DNA branch table from the web interface returned
HTTP 500, so the 16-leaf query is the reproducible public source currently
available from EvoNAPS.

## Run The Development Check

This command checks IQ-TREE wiring, target-branch selection and output parsing.
It is intentionally not described as a reconstruction of Fig. 2.

```bash
doc/satute-wiki/reconstruct_fig2_simulations.py \
  --iqtree ./build/iqtree3 \
  --outdir /tmp/iqtree-satute-paper-fig2-reconstruction \
  --reps 2 \
  --site-lengths 100,1000 \
  --branch-lengths 0.1,1.0,5.0
```

The development-check run committed into the wiki is:

```text
doc/satute-wiki/results/fig2_reduced_summary.tsv
doc/satute-wiki/results/fig2_reduced_detail.tsv
```

To launch the full paper-scale grid, use:

```bash
doc/satute-wiki/reconstruct_fig2_simulations.py \
  --iqtree ./build/iqtree3 \
  --outdir /tmp/iqtree-satute-paper-fig2-full \
  --reps 1000 \
  --paper-grid
```

The full grid is not run by default because it is a large simulation campaign.
A manuscript comparison requires the paper's Seq-Gen-based simulation details
and the EvoNAPS-derived branch-length sampling for the 16-taxon subtrees.
