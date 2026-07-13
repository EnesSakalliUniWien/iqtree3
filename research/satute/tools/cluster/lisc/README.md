# LiSC Full Head-To-Head Simulation

This folder prepares the paper-scale paired simulation for the University of
Vienna Life Science Compute Cluster. It uses Seq-Gen, IQ-TREE SatuTe and the
curated EvoNAPS branch-length table.

The default LiSC run is the paper-faithful Fig. 2 reproduction:

```text
simulation model: JC
evaluation model: JC
scenario set: fig2
tree cases: five_external,sixteen_internal
site lengths: 100,1000,10000
branch lengths: 0.1,0.2,0.3,0.4,0.5,0.8,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,7.5,10.0
replicates per point: 1000
formulas: dominant,eigenvalue_weighted (native IQ-TREE rows)
rows per replicate task: 8
```

This expands to 96,000 replicate tasks and 768,000 detail rows. The paper's
GTR analysis is a separate supplementary model-misspecification contract:
simulate under `GTR_PF06346` and evaluate under `GTR_PF06346`, `JC`, `K2P`, and
`F81` with `SCENARIO_SET=misspecification`.

Skewed GTR analyses are extension runs. Use `MODEL_PAIRS` for those so matched
simulation/evaluation designs stay explicit.

The default array has 1,000 shards, one per replicate index. Each paper-grid
shard therefore receives exactly one task from every combination of tree case,
site length and focal branch length: 96 tasks per shard. This avoids the 2:1
load imbalance produced by the earlier 768-shard layout.

## Setup On LiSC

From a network-enabled shell:

```bash
rsync -a --exclude .git --exclude build --exclude 'research/satute/artifacts/cluster/aws/full-run' \
  /Users/berksakalli/Projects/iq-tree-satute/ \
  lisc:/lisc/home/user/sakalli/projects/iq-tree-satute/
```

Then on LiSC:

```bash
ssh lisc
cd /lisc/home/user/sakalli/projects/iq-tree-satute
export WORK_DIR=/lisc/data/scratch/menche/sakalli/iq-tree-satute/full-head-to-head
bash research/satute/tools/cluster/lisc/build_iqtree_satute.sh
bash research/satute/tools/cluster/lisc/build_seqgen.sh
```

`build_seqgen.sh` uses `SEQGEN_SRC` when set. If the Seq-Gen source is already
available from the reanalysis project on LiSC, use:

```bash
SEQGEN_SRC=/lisc/home/user/sakalli/projects/satute-eigenvector-reanalysis/tools/seq-gen-1.3.4/source \
  bash research/satute/tools/cluster/lisc/build_seqgen.sh
```

## Fixed GTR20 protein reference

`fit_fixed_gtr20_reference.slurm` estimates the GTR20 exchangeabilities and
the `+G4` shape once on a large independent protein alignment and a supplied
reference topology. It writes checksums for every input and output. The fitted
model must then be copied into the experiment driver as a frozen `GTR20_REF`
alias; benchmark replicates must not invoke unconstrained `GTR20` fitting.

The fixed empirical protein matrix benchmark uses IQ-TREE AliSim so the same
backend supports LG, WAG, JTT, and Q.PFAM. Select it with `SIMULATOR=alisim`.
Run `run_protein_extension_canary.slurm` before submitting a larger protein
array; it checks all four matrices with and without fixed four-category gamma.

## Submit

Paper Fig. 2 JC-only run:

```bash
sbatch research/satute/tools/cluster/lisc/submit_full_head_to_head.slurm
```

Supplementary GTR model-misspecification run:

```bash
SIMULATION_MODELS=GTR_PF06346 \
EVALUATION_MODELS=GTR_PF06346,JC,K2P,F81 \
SCENARIO_SET=misspecification \
RUN_NAME=supp-misspec-gtr-pf06346-reps1000 \
  sbatch research/satute/tools/cluster/lisc/submit_full_head_to_head.slurm
```

Matched skewed-GTR extension run:

```bash
MODEL_PAIRS=GTR_SKEW_FREQ:GTR_SKEW_FREQ,GTR_SKEW_RATES:GTR_SKEW_RATES,GTR_SKEW_BOTH:GTR_SKEW_BOTH \
SCENARIO_SET=all \
RUN_NAME=extension-skewed-gtr-reps1000 \
  sbatch research/satute/tools/cluster/lisc/submit_full_head_to_head.slurm
```

Monitor:

```bash
squeue --me
sacct -j <jobid> --format=JobID,JobName,State,ExitCode,Elapsed,AllocCPUS,ReqMem
```

Each shard writes a readable `head_to_head_detail.tsv` and
`head_to_head_summary.tsv`. After a shard finishes successfully, the per-run
IQ-TREE and Seq-Gen files under `runs/` are compressed into:

```text
${OUT_ROOT}/run-archives/runs_shard_<index>_job_<jobid>.tar.gz
```

Shard 0 also writes `run_manifest.tsv` with the grid parameters and SHA-256
checksums for the exact IQ-TREE and Seq-Gen binaries, EvoNAPS table, experiment
driver and modular native SatuTe sources used by the run. The wrapper loads only
the pinned `GCCcore/14.3.0` runtime required by the IQ-TREE binary; it does not
load the SciPy/NumPy bundle.

The uncompressed `runs/` directory is removed only after the archive is written.

## Postprocess

After all array tasks complete:

```bash
export WORK_DIR=/lisc/data/scratch/menche/sakalli/iq-tree-satute/full-head-to-head
bash research/satute/tools/cluster/lisc/postprocess_full_head_to_head.sh
```

The default output root is:

```text
/lisc/data/scratch/menche/sakalli/iq-tree-satute/full-head-to-head/paper-fig2-jc-reps1000
```

Postprocessing combines shard TSVs, verifies the dominant statistic against
IQ-TREE SatuTe, checks the JC formula collapse, and generates the full-run
figures for every pair listed in `PLOT_MODEL_PAIRS`.
The dominant-statistic check uses `VERIFY_Z_TOLERANCE`, defaulting to `1e-4`,
because it compares an independent Python recomputation with IQ-TREE's C++
output.

For the supplementary GTR misspecification run, postprocess with:

```bash
RUN_NAME=supp-misspec-gtr-pf06346-reps1000 \
PLOT_MODEL_PAIRS=GTR_PF06346:GTR_PF06346,GTR_PF06346:JC,GTR_PF06346:K2P,GTR_PF06346:F81 \
PLOT_SCENARIO_SET=misspecification \
  bash research/satute/tools/cluster/lisc/postprocess_full_head_to_head.sh
```

For the matched skewed-GTR extension run, postprocess with:

```bash
RUN_NAME=extension-skewed-gtr-reps1000 \
PLOT_MODEL_PAIRS=GTR_SKEW_FREQ:GTR_SKEW_FREQ,GTR_SKEW_RATES:GTR_SKEW_RATES,GTR_SKEW_BOTH:GTR_SKEW_BOTH \
PLOT_SCENARIO_SET=all \
  bash research/satute/tools/cluster/lisc/postprocess_full_head_to_head.sh
```
