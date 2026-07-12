# Migration Inventory

The reorganization began with a dirty worktree and therefore preserves all
pre-existing user changes. No production file was reset or discarded.

## Pre-migration observations

- `doc/satute-wiki/results` occupied approximately 5.3 GB.
- Approximately 4.2 GB were ignored local runs and 1.1 GB were ignored cluster
  outputs.
- The root of `doc/satute-wiki` mixed more than thirty notes, drivers,
  renderers, build products, and manuscript files.
- `doc/satute-phase1` contained the native/reference regression suite but was
  ignored in its entirety.
- The manuscript consumed at least one figure from an ignored cluster path.

## Preservation decisions

- Native code remains at `tree/satute.cpp`.
- Raw local and cluster artifacts were moved, not deleted.
- Reviewed biological artifacts remain tracked under `artifacts/releases`.
- Historical documents were moved under `docs/archive` or `references`.
- Compatibility links preserve the old documentation and phase-1 paths.

## Resulting boundaries

- `src/` owns shared mathematical and artifact utilities.
- `experiments/` owns thin, numbered drivers and renderers.
- `tests/` owns exact, native/reference, and unit checks.
- `artifacts/local/` and `artifacts/cluster/` are ignored working data.
- `artifacts/releases/` contains compact reviewed evidence with checksums.
- `manuscript/` is self-contained and consumes only staged tracked figures.
- `references/` contains source papers, historical manuscript material, and
  external reference data.

## Post-migration gates

The experiment smoke profiles, exact checks, five-stage native regression
suite, structural audit, release checksum audit, high-resolution figure
renderer, and 17-page manuscript build all pass from the canonical workspace.
