# Experiment Architecture

The experiment suite uses one external seam:
`satute_analysis.experiments.ExperimentCatalog`. Maintenance commands discover,
validate, list, and run experiments through this interface. Manifest JSON is an
implementation detail: discovery converts it into immutable `Experiment`
values, validates all owned and project-relative paths, and returns defensive
configuration copies to workflows and release provenance.

Managed experiment directories supply workflow adapters at the internal seam.
An adapter owns profile-specific orchestration, reads its validated colocated
configuration, and delegates scientific work to that experiment's drivers. The
catalog loads only the adapter's
`run(context, *, profile, figures, study)` callable, so module-loading details
and command construction remain inside the implementation. Named studies make
distinct power and calibration designs visible without splitting one
scientific question across several top-level modules. Experiments that require
dataset-specific inputs declare
`execution: manual` and document direct driver commands instead of supplying a
workflow that cannot run.

## Ownership

- `experiments/<id>/` owns one scientific question, its manifest, workflow,
  drivers, renderers, and explanatory README.
- `src/satute_analysis/` owns reusable mathematics, schemas, parsing,
  statistics, artifact handling, and catalog discovery.
- `tests/` verifies public interfaces and observable scientific behavior.
- `tools/maintenance/` contains thin command-line adapters for the catalog.
- `tools/cluster/` contains environment-specific execution adapters.
- `artifacts/local/` and `artifacts/cluster/` contain disposable run outputs.
- `artifacts/releases/` contains reviewed, checksummed evidence.
- `manuscript/` consumes promoted release evidence, never live run directories.

## Dependency direction

```text
./satute
        |
        v
ExperimentCatalog
        |
        v
named study ---> experiment workflow.py ---> experiment driver / renderer
        |                              |
        +------------------------------+
                       |
                       v
             satute_analysis shared modules
```

Scientific drivers do not import peer experiment scripts. Code reused by two
experiments belongs in `satute_analysis`; environment-specific shell behavior
belongs in `tools`.

## Interface invariants

- Experiment identifiers match their directory names.
- Driver and managed-workflow paths remain inside the owning experiment module.
- Every numbered experiment has a colocated manifest and README.
- Study identifiers are unique and every default study exists.
- Output paths are dedicated children of `artifacts/local/` or
  `artifacts/cluster/`; release paths are dedicated children of
  `artifacts/releases/`.
- Output directories are unique across discovered experiments.
- Consumers are unique non-empty names.
- Configuration schema versions are positive integers.

Callers and tests verify these invariants through `ExperimentCatalog`; they do
not parse manifests independently.

## Adding an experiment

1. Create `experiments/NNN_descriptive_name/`.
2. Add a valid `experiment.json` with at least one named study, a driver,
   README, and either a managed `workflow.py` or `execution: manual`.
3. Keep the driver CLI independently usable with `--help`.
4. Add behavioral tests under the matching `tests/` category.
5. Run `./satute check audit`; discovery is automatic and no central registry edit is
   required.
