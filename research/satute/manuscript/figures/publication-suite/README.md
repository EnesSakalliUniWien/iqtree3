# Publication-suite figure index

These tracked PDF panels are copied from the verified LiSC outputs and are the
small main-text subset.  The complete figure trees remain in the ignored
cluster artifact directories, with one independent output directory for every
model pair and decision rule.

- Source postprocessor: immutable LiSC checkout `7b264fe9c5ef`.
- Simulation binary SHA: `a5982477099d6dc64e9c0182a37255048b7f6e28c810d6dc861c4c1a6fbd2646`.
- Publication curves: 1,000 replicates per cell, 100/250/1,000 sites,
  branch lengths 4--12, five-external and sixteen-internal trees.
- Each R output directory contains independent dominant, eigenvalue-weighted,
  and paired-difference 2D PDF/SVG/450-dpi PNG plus vector 3D PDF and
  self-contained HTML.
- Calibration figures use the exact-null (5,000 replicates, 2,500 training)
  and mixed-null (2,000 replicates) truth-labelled designs and report
  unadjusted, taxon-Bonferroni, and BY rules.
- The tracked calibration panels are copied from
  `fdr-calibration-reps5000-2000-corrected-ties`. Its training-derived
  thresholds include complete adjusted-p-value tie blocks only, matching the
  implemented `p <= threshold` decision rule.
- The tracked protein panels show the taxon-Bonferroni paired differences for
  fixed LG+G4, WAG+G4, JTT+G4, and Q.pfam+G4. The complete protein artifact
  tree contains every model and all three decision rules in 2D and 3D formats.
