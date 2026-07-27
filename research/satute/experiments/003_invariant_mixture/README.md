# 003 - Invariant-mixture correction

Evaluates the counterfactual and corrected `+I` likelihood constructions on
paired exact-pattern multinomial alignments.

- Driver: `run.py`
- Full renderer: `render_full.py`
- Direct-shift renderer: `render_shift.py`
- Manifest and run configuration: `experiment.json`
- Common workflow adapter: `workflow.py`
- Reviewed release: `../../artifacts/releases/invariant-mixture-v1/`

Smoke test:

```bash
../../satute data --profile smoke
```

Full run:

```bash
../../satute experiment 003_invariant_mixture --profile full
../../satute experiment 003_invariant_mixture --profile full --figures
```
