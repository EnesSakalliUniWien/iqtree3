# Enhanced Sliding-Window Biological Rerun

This directory contains the curated nucleotide Tree-of-Life sliding-window
rerun from 2026-06-19. It covers the manuscript 16S rRNA `rRNA_based_3D_tree`
case under the fitted GTR+F+G4 model. The focal branches are the internal
branch leading to Eukaryota and the external branch leading to
`Saccharomyces_cerevisiae`.

The curated comparison keeps only the published SatuTe statistic and the
eigenvalue-weighted extension. Method-specific outputs are stored under:

```text
original_satute/rRNA_based_3D_tree/
sliding_window/dominant/rRNA_based_3D_tree/
sliding_window/eigenvalue_weighted/rRNA_based_3D_tree/
components/dominant/rRNA_based_3D_tree/
components/eigenvalue_weighted/rRNA_based_3D_tree/
```

The summary table used by the manuscript is:

```text
comparison_to_manuscript/original_vs_enhanced_summary.tsv
```

The main comparison plots are:

```text
comparison_to_manuscript/original_vs_enhanced_sliding_window_curves.pdf
comparison_to_manuscript/original_vs_enhanced_saturated_windows.pdf
```

The protein LG+G4 Tree-of-Life case was not recomputed with the enhanced
statistic in this run.
