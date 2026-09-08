# SatuTe Research Workspace

This directory contains the research, validation, experiments, release
artifacts, and manuscript supporting the SatuTe integration in IQ-TREE.
Production code remains in [`../../tree/satute.cpp`](../../tree/satute.cpp).

## Start Here

- [`PROJECT_STATUS.md`](PROJECT_STATUS.md): current implementation and evidence.
- [`experiments/`](experiments/): self-contained experiment modules and their manifests.
- [`docs/overview/experiment-architecture.md`](docs/overview/experiment-architecture.md):
  module ownership and dependency direction.
- [`artifacts/releases/`](artifacts/releases/): reviewed, compact evidence.
- [`manuscript/main.tex`](manuscript/main.tex): modular manuscript entry point.
- [`tests/`](tests/): exact, native, and regression validation.
- [`docs/`](docs/): theory, experiment notes, paper map, and archive.

## Commands

Run from this directory through the sole workspace interface:

```bash
./satute bootstrap
./satute doctor
./satute status
./satute data --profile smoke
./satute experiment 003_invariant_mixture
./satute experiment 003_invariant_mixture --profile full --figures
./satute cluster plan
./satute manuscript
./satute check audit
```

`./satute data --profile smoke` runs compact native-backed profiles for
experiments 001--004. The launcher always uses the single project environment
at `.venv`; internal Python scripts and Make targets are implementation
adapters, not additional supported interfaces.
Experiment 005 intentionally requires explicit biological alignment inputs.
Use `PROFILE=full` only for a deliberate full rerun.

The versioned expanded recalculation contract is
`campaigns/2026_expanded_recalculation.json`. Inspect it with
`./satute cluster plan`; submit its cluster-owned runs from a clean,
source-matched LiSC checkout with `./satute cluster submit`. The plan preserves
the original JC grid as a named baseline and appends longer branches without
changing the reconstructed original-paper comparison.

Experiment runs write to `artifacts/local/` by default and should leave Git
clean. Reviewed summaries and final figures are promoted explicitly to
`artifacts/releases/<experiment>/<version>/`.

## Workspace Boundaries

| Path | Purpose | Git policy |
|---|---|---|
| `src/` | Shared Python analysis library | tracked |
| `tests/` | Exact, native, and regression tests | tracked |
| `experiments/` | Self-contained manifests, workflows, drivers, and renderers | tracked |
| `artifacts/releases/` | Reviewed summaries and final figures | tracked |
| `artifacts/local/` | Workstation runs | ignored |
| `artifacts/cluster/` | LiSC/AWS raw outputs and shards | ignored |
| `manuscript/` | Manuscript source and staged assets | tracked except build |
| `references/` | Papers and historical source material | tracked |
| `docs/archive/` | Superseded plans and compatibility notes | tracked |

Experiments are discovered from their colocated `experiment.json` manifests;
the common execution interface is implemented by
`src/satute_analysis/experiments.py`.
The canonical manuscript order is defined by `manuscript/main.tex`; Results
and Online Methods are further split into focused modules under
`manuscript/sections/results/` and `manuscript/sections/methods/`.

## Command Reference

Every subcommand below is `./satute <command> ...`, run from this directory.

### `bootstrap`

Creates/updates the sole `.venv` environment (editable install of
`satute-analysis` and pinned dependencies), then runs `doctor`. No options.

### `doctor [--json]`

Verifies the entry point, interpreter, virtual environment, and IQ-TREE
binary are consistent (fails if more than one venv-like directory exists in
the repository, or the running interpreter isn't `.venv`'s). `--json` prints
a machine-readable report instead of the text summary. Exit code is non-zero
when the environment is inconsistent.

### `status`

Reports current experiment and artifact status. No options.

### `data [--profile {smoke,full}] [--map] [--only EXPERIMENT/STUDY] [--plan] [--skip-rate-category-check] [--figures] [--figure-mode {render,index-existing}]`

Plans or runs managed data generation across every registered experiment
study (excluding `scale: manual` studies, e.g. experiment 005's biological
studies, which require explicit alignment inputs).

- `--profile {smoke,full}` (default `smoke`): compact wiring-validation
  replicate counts vs. the full declared evidence grid.
- `--map`: print the complete simulation/analysis/artifact-lineage map
  instead of running under `--profile`.
- `--only EXPERIMENT/STUDY` (repeatable): restrict to specific studies, e.g.
  `--only 002_relative_weighting/gtr_power`.
- `--plan`: print the run plan without executing any driver.
- `--skip-rate-category-check`: skip the native rate-category cross-check.
- `--figures`: after data generation, also build the canonical publication
  power figures from the resulting case map.
- `--figure-mode {render,index-existing}` (default `render`, only used with
  `--figures`): render figures from source tables, or just index
  already-approved manuscript PDFs without re-rendering.

### `experiment EXPERIMENT [--study STUDY] [--profile {smoke,full}] [--figures]`

Runs one registered experiment's study directly (`EXPERIMENT` is the
directory name, e.g. `002_relative_weighting`).

- `--study STUDY`: study id from that experiment's `experiment.json`
  (defaults to the manifest's `default_study`).
- `--profile {smoke,full}` (default `smoke`).
- `--figures`: **renders only** — reads the study's already-written
  `head_to_head_summary.tsv` and produces figures from it; it does **not**
  run data generation itself. Run the same command once *without*
  `--figures` first, then rerun *with* `--figures` to render from that
  output. Passing `--figures` before any prior run fails with
  `FileNotFoundError` on the missing summary table.

```bash
./satute experiment 002_relative_weighting --study gtr_power --profile smoke
./satute experiment 002_relative_weighting --study gtr_power --profile smoke --figures
```

### `publication-figures [--plan] [--index-existing]`

Builds the canonical publication figure suite under `artifacts/publication/`
from the full case map.

- `--plan`: print the plan without rendering.
- `--index-existing`: index already-rendered/approved figures instead of
  rendering from source tables.

### `cluster {plan,submit} [--campaign PATH] [--only RUN_ID]`

Plans or submits the versioned LiSC recalculation campaign.

- `--campaign PATH` (default `campaigns/2026_expanded_recalculation.json`).
- `--only RUN_ID` (repeatable): limit to specific campaign run ID(s).
- `plan` prints the resolved plan; `submit` actually launches it and should
  only be run from a clean, source-matched LiSC checkout.

### `check [{all,power,figures,publication-figures,rate-categories,audit}]`

Runs a reproducible validation suite (default `all`).

- `all`: full unit-test suite, exact-weighted-likelihood test, a
  2000-replicate exact-null check, the native phase-1 checks, and the repo
  hygiene audit.
- `power`: the unit tests from `test_cluster_campaign.py` onward (skipping
  the experiment-catalog and workspace-interface tests), the
  exact-weighted-likelihood test, and the native rate-category cross-check.
- `figures`: figure/manuscript-integration checks (`audit_power_plots.py`,
  `verify_figure_integration.py`, `journal_artwork_gate.py`,
  `verify_equation_integration.py`).
- `publication-figures`: publication-figure contract tests plus the
  power-plot audit.
- `rate-categories`: native rate-category cross-check only.
- `audit`: repo hygiene audit only.

### `manuscript`

Compiles `manuscript/main.tex` and `manuscript/supplementary.tex` with
`latexmk`, then copies the results to `manuscript/manuscript.pdf` and
`manuscript/supplementary-information.pdf`. No options; requires a working
LaTeX toolchain.

### `workflow-figure`

Renders the LaTeX workflow figure (Figure 1). No options.

### `release EXPERIMENT VERSION`

Creates a reviewed experiment release manifest under
`artifacts/releases/<experiment>/<version>/`. Both arguments are required
positionals.

### `developer ADAPTER [-- ARGS...]`

Runs one allow-listed low-level script directly through the canonical
environment, bypassing the experiment/study abstraction. Anything after
`ADAPTER` (optionally preceded by `--`) is forwarded verbatim to the
underlying script:

```bash
./satute developer paper-head-to-head -- --branch-lengths 0.5,1.0,2.0,3.0 --reps 20 \
  --site-lengths 1000 --tree-cases five_external \
  --model-pairs GTR_PF06346:GTR_PF06346 --scenario-set fig2 --simulator alisim \
  --iqtree ../../build/iqtree3 --outdir artifacts/local/demo
```

| adapter | script |
|---|---|
| `biological-analysis` | `experiments/005_biological_windows/analyze_dataset.py` |
| `formula-comparison` | `tests/native/compare_formulas.py` |
| `invariant-full-figure` | `experiments/003_invariant_mixture/render_full.py` |
| `invariant-shift-figure` | `experiments/003_invariant_mixture/render_shift.py` |
| `invariant-simulation` | `experiments/003_invariant_mixture/run.py` |
| `paper-head-to-head` | `experiments/002_relative_weighting/run.py` |
| `reference` | `tests/native/satute_reference.py` |
| `reproduce-scientific-results` | `tools/maintenance/reproduce_scientific_results.py` |
| `saturation-null-check` | `tests/exact/test_saturation_null.py` |
| `saturation-simulation` | `experiments/004_saturation_diagnostics/run.py` |
| `validate-weighted-formula` | `tests/native/validate_weighted_formula.py` |

### `clean [--dry-run] [--report] [--purge ENTRY ...]`

With no options, removes only output that one documented command rebuilds in
full: `manuscript/build/`, `artifacts/local/smoke/`, `artifacts/local/check/`,
`.pytest_cache/`, and any `__pycache__` outside `.venv`. It reports what each
removal costs to rebuild, and the total reclaimed.

- `--dry-run`: list what would be removed, with sizes, and delete nothing.
- `--report`: inventory every entry under `artifacts/local/` with its size and
  provenance — regenerable, declared by an experiment manifest, or ad hoc.
- `--purge ENTRY [...]`: delete named entries under `artifacts/local/`.

`clean` deliberately never touches ad hoc entries. Most of what accumulates in
`artifacts/local/` is one-off exploratory runs that no registered command
recreates, so deleting them loses the results permanently — they have to be
named explicitly via `--purge`. Paths outside `artifacts/local/`, and
`artifacts/local` itself, are refused.

```bash
./satute clean --report
./satute clean --purge skewed-gtr-sites-reps500
```
