# Experiment Modules

Each numbered directory is one self-contained experiment module. Its owned
files follow the same small interface:

| File | Ownership |
|---|---|
| `experiment.json` | Identity, question, status, consumers, paths, and run configuration |
| `workflow.py` | Smoke/full orchestration behind `make experiment` and `make figures` |
| `run.py` or the declared driver | Scientific computation and its command-line interface |
| `render*.py` / `render*.R` | Figure or table rendering for that experiment |
| `README.md` | Scientific intent, evidence contract, and direct-use notes |

The directories are intentionally numbered because execution order and
historical dependency matter. They are discovered from their colocated
manifests; there is no second central registry to keep synchronized.

Shared mathematics, schemas, artifact handling, and experiment discovery live
under `src/satute_analysis/`. Cluster and maintenance adapters live under
`tools/`; generated outputs never live beside experiment source.
Smoke outputs are isolated under
`artifacts/local/smoke/<experiment-id>/schema-<config-schema>` so an old result
cannot be resumed under a newer table contract.

Use the common interface from the research workspace:

```bash
make experiment EXP=003_invariant_mixture PROFILE=smoke
make figures EXP=003_invariant_mixture PROFILE=full
```

Experiment 005 requires explicit biological inputs and therefore documents
direct driver commands instead of inventing defaults.
