# Experiment Architecture

The experiment suite uses one external seam:
`satute_analysis.experiments.ExperimentCatalog`. Maintenance commands discover,
validate, list, and run experiments through this interface. Manifest JSON is an
implementation detail: discovery converts it into immutable `Experiment`
values, validates all owned and project-relative paths, and returns defensive
configuration copies only for release provenance.

Each numbered experiment directory supplies one workflow adapter at the
internal seam. The adapter owns profile-specific orchestration and delegates
scientific work to that experiment's drivers. The catalog loads only the
adapter's `run(context, *, profile, figures)` callable, so module-loading details
and command construction remain inside the implementation. This gives five
real adapters behind one interface.

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
Make / maintenance CLI
        |
        v
ExperimentCatalog
        |
        v
experiment workflow.py  --->  experiment driver / renderer
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
- Driver and workflow paths remain inside the owning experiment module.
- Output and release paths remain inside the SatuTe research project.
- Output directories are unique across discovered experiments.
- Consumers are unique non-empty names.
- Configuration schema versions are positive integers.

Callers and tests verify these invariants through `ExperimentCatalog`; they do
not parse manifests independently.

## Adding an experiment

1. Create `experiments/NNN_descriptive_name/`.
2. Add a valid `experiment.json`, `workflow.py`, driver, and README.
3. Keep the driver CLI independently usable with `--help`.
4. Add behavioral tests under the matching `tests/` category.
5. Run `make audit`; discovery is automatic and no central registry edit is
   required.
