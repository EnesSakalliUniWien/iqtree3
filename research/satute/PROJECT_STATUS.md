# Project Status

## Production

- Native implementation: `../../tree/satute.cpp`.
- Native formulas: `dominant` and `eigenvalue_weighted`.
- Pooled branch p-values use separate Benjamini-Yekutieli FDR families for
  each native formula.
- Monotonic model ruler: `InformationFraction` and `SaturationIndex`.

## Current Mathematical Result

For an explicit invariant category, the infinite-branch null uses the
invariant joint likelihood. The invariant component contributes no
finite-versus-infinite-branch numerator signal. Exact enumeration verifies the
refined likelihood identity and null centering.

## Current Empirical Release

`artifacts/releases/invariant-mixture-v1/` contains the reviewed dense paired
`+I+G4`/`+I+R4` grid:

- 27 branch lengths;
- 5 alignment lengths;
- 20,000 null and 20,000 alternative replicates per cell;
- 10.8 million paired simulated alignments;
- maximum refined power gain: 0.16315;
- refined mean nominal null rejection: 0.0492;
- counterfactual mean nominal null rejection: 0.9002.

## Validation

- Exact likelihood and null-centering checks: passing.
- Native/reference comparisons: passing within `2e-4`.
- Five-stage native regression suite: passing.
- Registry driver CLI, cluster shell, compatibility-link, and release-checksum
  audit: passing.
- Compact experiment smoke profiles 001--004: passing.
- Manuscript build: passing.

Run `make status` for filesystem and Git state and `make check` for validation.
