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
4. [04-new-eigenvalue-statistic.md](04-new-eigenvalue-statistic.md) separates
   the published dominant-mode statistic from the new all-mode and
   eigenvalue-decay weighted statistics.
5. [head_to_head_satute_simulations.py](head_to_head_satute_simulations.py)
   is the paired simulation driver: one simulated alignment is evaluated by
   the dominant, all-mode and eigenvalue-weighted tests under the same
   scenario and threshold.
6. [short-addendum-paper.tex](short-addendum-paper.tex) is a manuscript draft
   for the eigenvalue-weighted extension. It uses the SatuTe notation and does
   not present the reduced development checks as empirical evidence. The
   compiled PDF is [short-addendum-paper.pdf](short-addendum-paper.pdf).
7. [five-page-explainer.md](five-page-explainer.md) is the longer working
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
  --models JC \
  --simulator alisim
```

The output has one row per replicate, analysis scenario and formula. The same
alignment and target branch are evaluated by the dominant, all-mode and
eigenvalue-weighted statistics.

The verified development-check outputs are:

```text
doc/satute-wiki/results/head_to_head_development_detail.tsv
doc/satute-wiki/results/head_to_head_development_summary.tsv
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
  --models JC
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
