# Independent Enhanced Biological Sliding-Window Rerun

This directory contains the curated subset of the 2026-06-19 independent
eigenvalue-weighted sliding-window rerun. It keeps the manuscript SatuTe rerun,
the native dominant calculation and the eigenvalue-weighted calculation. The
full component CSVs are not copied here; the curated tables retain provenance
paths to the source run.

## Datasets

The protein example uses the manuscript `protein_based_2D_tree` Tree-of-Life
case:

```text
trees/protein_based_2D_ToL.treefile
```

The nucleotide example uses the manuscript `rRNA_based_3D_tree` Tree-of-Life
case:

```text
trees/rRNA_based_3D_ToL.treefile
```

Both analyses use 36-site sliding windows and the same two target branches:
the branch leading to Eukaryota and the branch leading to yeast.

## Window Counts

The focused published SatuTe versus eigenvalue-weighted comparison is:

```text
comparison_to_manuscript/published_vs_eigenvalue_weighted_summary.tsv
```

In the protein 2D tree, the branch leading to Eukaryota has 2,561 windows; the
published SatuTe rerun marks 138 windows as saturated and the
eigenvalue-weighted calculation marks 5. The yeast branch also has 2,561
windows; the corresponding counts are 54 and 8.

In the 16S rRNA 3D tree, the branch leading to Eukaryota has 1,912 windows; the
published SatuTe rerun marks 822 windows as saturated and the
eigenvalue-weighted calculation marks 699. The yeast branch also has 1,912
windows; the corresponding counts are 70 and 29.

## Figures

The manuscript-ready comparison figures are:

```text
comparison_to_manuscript/published_vs_eigenvalue_weighted_curves.pdf
comparison_to_manuscript/published_vs_eigenvalue_weighted_saturated_windows.pdf
```

They can be regenerated with:

```bash
python3 doc/satute-wiki/plot_independent_biological_comparison.py
```
