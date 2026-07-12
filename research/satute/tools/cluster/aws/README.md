# AWS Full-Run Notes

This folder is for running the paper-scale, paired SatuTe simulation on AWS.
It does not use development-check data.

## Required Inputs

The full paper-faithful run requires:

1. A Linux build of `iqtree3` from this repository.
2. Seq-Gen v1.3.4 or compatible `seq-gen` executable.
3. An EvoNAPS branch-length table with columns:
   - `type`, `kind`, or `branch_type`
   - `length`, `branch_length`, `blen`, or `bl`

The table must contain internal and external branch lengths from the EvoNAPS
empirical pools, filtered by the script to `0.04 <= length <= 0.2`. The
official EvoNAPS web download uses `BRANCH_TYPE` with `i/e` values and `BL` for
length; this format is accepted directly.

The recovered public branch-length source for the 16-taxon protocol is:

```text
research/satute/references/evonaps/evonaps_16taxon_branch_lengths.tsv
```

It was downloaded from the EvoNAPS tree search result page for DNA trees with
16 leaves and maximum internal/external branch lengths no greater than `0.2`.
After the manuscript interval filter, it contains 92 internal and 233 external
branch lengths.

## Per-Shard Command

Each shard runs a disjoint subset of replicate tasks:

```bash
ROOT_DIR=/work/iq-tree-satute \
PYTHON_BIN=/usr/bin/python3 \
IQTREE_BIN=/work/iq-tree-satute/build/iqtree3 \
SEQGEN_BIN=/work/bin/seq-gen \
EVONAPS_BRANCH_LENGTHS=/work/input/evonaps_16taxon_branch_lengths.tsv \
OUTDIR=/work/output/satute-head-to-head \
REPS=1000 \
SHARD_INDEX=${AWS_BATCH_JOB_ARRAY_INDEX} \
SHARD_COUNT=200 \
SIMULATION_MODELS=JC \
EVALUATION_MODELS=JC \
SCENARIO_SET=fig2 \
TREE_CASES=five_external,sixteen_internal \
  /work/iq-tree-satute/research/satute/tools/cluster/aws/run_full_head_to_head_shard.sh
```

This is the paper-faithful Fig. 2 setting. Fig. 2 is JC-only. The paper's GTR
analysis is the separate Supplementary Fig. model-misspecification experiment
using the EvoNAPS PF06346 model alias `GTR_PF06346`; do not add it to the Fig.
2 rerun.

For the supplementary GTR model-misspecification contract, use:

```bash
SIMULATION_MODELS=GTR_PF06346 \
EVALUATION_MODELS=GTR_PF06346,JC,K2P,F81 \
SCENARIO_SET=misspecification \
TREE_CASES=five_external,sixteen_internal \
  /work/iq-tree-satute/research/satute/tools/cluster/aws/run_full_head_to_head_shard.sh
```

For matched skewed-GTR extension runs, use exact model pairs:

```bash
MODEL_PAIRS=GTR_SKEW_FREQ:GTR_SKEW_FREQ,GTR_SKEW_RATES:GTR_SKEW_RATES,GTR_SKEW_BOTH:GTR_SKEW_BOTH \
SCENARIO_SET=all \
TREE_CASES=five_external,sixteen_internal \
  /work/iq-tree-satute/research/satute/tools/cluster/aws/run_full_head_to_head_shard.sh
```

If a shard is interrupted after writing partial output, rerun the same shard
with `RESUME=1`. The driver appends to the existing detail TSV and skips
replicate tasks that already have the expected row count for the selected
scenario set and model pair.

## Combine, Verify And Plot

After all shards finish and their output folders are available on one machine:

```bash
ROOT_DIR=/work/iq-tree-satute \
PYTHON_BIN=/usr/bin/python3 \
SHARD_ROOT=/work/output/satute-head-to-head \
OUTDIR=/work/output/full-run-results \
FIGURE_DIR=/work/output/full-run-figures \
EXPECTED_REPS=1000 \
TREE_CASES=five_external,sixteen_internal \
PLOT_MODEL_PAIRS=JC:JC \
PLOT_SCENARIO_SET=fig2 \
  /work/iq-tree-satute/research/satute/tools/cluster/aws/combine_verify_plot_full.sh
```

This command:

1. Combines all `head_to_head_detail.tsv` shard files.
2. Recomputes `head_to_head_summary.tsv`.
3. Verifies that the dominant formula matches IQ-TREE SatuTe and that JC
   collapses across formulas.
4. Generates model-specific manuscript plots:
   - `figure_simulated_data_main_reproduction_sim_JC__eval_JC.svg/pdf`
   - `figure_head_to_head_formula_comparison_sim_JC__eval_JC.svg/pdf`

The plotter refuses incomplete paper-scale input unless `--allow-incomplete` is
passed manually. Do not use `--allow-incomplete` for manuscript figures.

For long-running EC2 jobs, `watch_combine_full.sh` can be launched under
`nohup`; it waits for all shard workers to finish and then runs the combine,
verification and plotting step.
