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

The former `doc/satute-wiki` and `doc/satute-phase1` paths are compatibility
links for directory-level navigation and should not be used for new commands.
Experiments are discovered from their colocated `experiment.json` manifests;
the common execution interface is implemented by
`src/satute_analysis/experiments.py`.
The canonical manuscript order is defined by `manuscript/main.tex`; Results
and Online Methods are further split into focused modules under
`manuscript/sections/results/` and `manuscript/sections/methods/`.
