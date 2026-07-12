# 003 - Invariant-mixture correction

Evaluates the counterfactual and corrected `+I` likelihood constructions on
paired exact-pattern multinomial alignments.

- Driver: `run.py`
- Full renderer: `render_full.py`
- Direct-shift renderer: `render_shift.py`
- Config: `../../config/experiments/003_invariant_mixture.yaml`
- Reviewed release: `../../artifacts/releases/invariant-mixture-v1/`

Smoke test:

```bash
make -C ../.. smoke
```

Full run:

```bash
make -C ../.. experiment EXP=003_invariant_mixture PROFILE=full
make -C ../.. figures EXP=003_invariant_mixture
```
