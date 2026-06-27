# Paper Simulation Target And Smoke Test

The development script is:

```text
doc/satute-wiki/reconstruct_fig2_simulations.py
```

The paired head-to-head script for comparing formulas is:

```text
doc/satute-wiki/head_to_head_satute_simulations.py
```

It implements the four Fig. 2 analysis scenarios:

```text
true_tree_fixed_lengths
true_topology_ml_lengths
ml_tree_unadjusted
ml_tree_bonferroni
```

The five-taxon manuscript tree tests the `B | A1,A2,A3,A4` split. Taxon `B` is
connected by the focal branch `AB` to the root `A` of a balanced four-taxon
subtree whose internal and external branches are fixed at `0.2`:

```text
(((A1:0.2,A2:0.2):0.2,(A3:0.2,A4:0.2):0.2):0.0,B:L);
```

The 16-taxon tree is a balanced split:

```text
A1,...,A8 | B1,...,B8
```

The focal branch is the central branch between those two subtrees.

The reduced outputs committed in this folder are development checks. They use
small replicate counts and do not reproduce the paper simulations.

The manuscript archive in `sources/SatuTe_manuscript_2024_08` gives the
pre-publication Nature-format protocol. It clarifies that the five-taxon
simulation is the `B | A1,A2,A3,A4` split: taxon `B` is connected by branch
`AB` to the root `A` of a balanced four-taxon subtree. It also states that
Seq-Gen v1.3.4 was used, and that the 16-taxon subtree branch lengths were
sampled from EvoNAPS internal and external branch-length distributions.

The public EvoNAPS web interface exposes a branch-length download for tree
search results. The recovered 16-taxon branch source is:

```text
doc/satute-wiki/sources/evonaps/evonaps_16taxon_branch_lengths.tsv
```

The file was downloaded from the official EvoNAPS tree search result page for
DNA trees with 16 leaves and maximum internal/external branch length no greater
than `0.2`. The driver then applies the manuscript interval filter
`0.04 <= length <= 0.2`, leaving 92 internal and 233 external candidate lengths.
An unfiltered all-DNA branch download from the EvoNAPS web interface returned
HTTP 500, so this 16-leaf query is the reproducible public source currently
available.

## Smoke-Test Run

The development-check run used:

```text
reps = 2
site lengths = 100, 1000
branch lengths = 0.1, 1.0, 5.0
```

The output files are:

```text
doc/satute-wiki/results/fig2_reduced_detail.tsv
doc/satute-wiki/results/fig2_reduced_summary.tsv
```

These files only verify command-line wiring, target branch identification,
output parsing and plotting. They are too small, and the simulator/tree details
are too different, to support claims about the published power curves.

## Head-To-Head Formula Comparison

The enhancement is evaluated on paired replicates. For each simulated
alignment and each analysis scenario, the same target branch is evaluated by
the published dominant SatuTe statistic and the eigenvalue-weighted statistic.
The same unadjusted or Bonferroni-adjusted threshold is then applied to both
formulas. This produces paired decision rows and discordance counts between the
published and enhanced tests.

Smoke-test command:

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

Paper-faithful command shape:

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

For full EC2 execution, use `doc/satute-wiki/aws/run_full_head_to_head_shard.sh`
with `SEQGEN_BIN` and, for the 16-taxon case, `EVONAPS_BRANCH_LENGTHS` pointing
to the recovered EvoNAPS TSV.

This is the Fig. 2 simulation design. The main paper simulations are
JC-only: Seq-Gen v1.3.4 under JC, IQ-TREE/SatuTe under JC, 1,000 replicates per
point, site lengths `100,1000,10000`, the full 16-value branch grid, and the
four analysis scenarios.

The GTR experiment in the paper is not Fig. 2. It is the Supplementary Fig.
model-misspecification experiment: simulate under the EvoNAPS PF06346 GTR model
(`GTR_PF06346`) with rates `0.6676,3.7807,4.2833,0.5354,0.8718,1.0` and
frequencies `0.125,0.436,0.191,0.245`, then evaluate under correctly specified
GTR and under misspecified JC, K2P and F81. An eigenvector enhancement should
be evaluated on these paired designs rather than added as a placeholder GTR
model in the Fig. 2 reproduction.

The corresponding supplementary command shape is:

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

Skewed GTR models are useful for stress-testing the eigenvalue-weighted
calculation, but they are not part of the published simulation reproduction.
Use exact model pairs for this extension:

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

## Full Grid Command

The full paper-scale command is:

```bash
doc/satute-wiki/reconstruct_fig2_simulations.py \
  --iqtree ./build/iqtree3 \
  --outdir /tmp/iqtree-satute-paper-fig2-full \
  --reps 1000 \
  --paper-grid
```

That command is not run by default. It expands to two tree cases,
three sequence lengths, sixteen focal branch lengths, one thousand replicates,
and four analysis scenarios per replicate. Before using it for manuscript
claims, prefer the paired head-to-head driver above, because it uses Seq-Gen
and can consume the recovered EvoNAPS branch-length TSV.
