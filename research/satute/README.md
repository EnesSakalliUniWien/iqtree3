# SatuTe Research Workspace

This directory contains the research, validation, experiments, release
artifacts, and manuscript supporting the SatuTe integration in IQ-TREE.
Production code remains in [`../../tree/satute.cpp`](../../tree/satute.cpp).

## Start Here

- [`PROJECT_STATUS.md`](PROJECT_STATUS.md): current implementation and evidence.
- [`experiments/registry.yaml`](experiments/registry.yaml): authoritative experiment registry.
- [`artifacts/releases/`](artifacts/releases/): reviewed, compact evidence.
- [`manuscript/main.tex`](manuscript/main.tex): modular manuscript entry point.
- [`tests/`](tests/): exact, native, and regression validation.
- [`docs/`](docs/): theory, experiment notes, paper map, and archive.

## Commands

Run from this directory:

```bash
make status
make smoke
make check
make experiment EXP=003_invariant_mixture
make figures EXP=003_invariant_mixture
make manuscript
make audit
```

`make smoke` runs compact native-backed profiles for experiments 001--004.
Experiment 005 intentionally requires explicit biological alignment inputs.
Use `PROFILE=full` only for a deliberate full rerun.

Experiment runs write to `artifacts/local/` by default and should leave Git
clean. Reviewed summaries and final figures are promoted explicitly to
`artifacts/releases/<experiment>/<version>/`.

## Workspace Boundaries

| Path | Purpose | Git policy |
|---|---|---|
| `src/` | Shared Python analysis library | tracked |
| `tests/` | Exact, native, and regression tests | tracked |
| `experiments/` | Thin experiment drivers and configs | tracked |
| `artifacts/releases/` | Reviewed summaries and final figures | tracked |
| `artifacts/local/` | Workstation runs | ignored |
| `artifacts/cluster/` | LiSC/AWS raw outputs and shards | ignored |
| `manuscript/` | Manuscript source and staged assets | tracked except build |
| `references/` | Papers and historical source material | tracked |
| `docs/archive/` | Superseded plans and compatibility notes | tracked |

The former `doc/satute-wiki` and `doc/satute-phase1` paths are compatibility
links for directory-level navigation and should not be used for new commands.
The canonical manuscript order is defined by `manuscript/main.tex`; Results
and Online Methods are further split into focused modules under
`manuscript/sections/results/` and `manuscript/sections/methods/`.
