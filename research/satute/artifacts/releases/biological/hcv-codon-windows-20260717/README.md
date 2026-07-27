# HCV codon-position sliding-window example

This release contains the second empirical SatuTe example used in the
manuscript. It is an illustrative reference-genome panel, not an epidemiological
sample or an analysis of within-genotype diversity.

## Reference panel

Ten complete hepatitis C virus (HCV) records were retrieved from NCBI
Nucleotide by accession version: `NC_004102.1`, `AJ238799.1`, `AB047639.1`,
`D10988.1`, `D17763.1`, `D49374.1`, `Y11604.1`, `AF064490.1`, `Y12083.1` and
`EF108306.2`. The combined GenBank flat file is retained in `source/`.

The primary polyprotein CDS was translated, aligned with MAFFT 7.526 and
back-translated to a codon-preserving nucleotide alignment. H77
(`NC_004102.1`) mature-peptide annotations define Core--NS5B coordinates. A
full-coding topology was inferred under GTR+F+G4; the subtype pairs for
genotypes 1, 2 and 3 define the three tested internal branches.

## Window design and validation

Codon positions 1+2 and position 3 were concatenated separately. Every
36-site window, advanced by 12 sites, was refitted under GTR+F+G4 on the fixed
topology. Thus positions 1+2 and position 3 have the same statistical sample
size but span 18 and 36 aligned codons, respectively. Saturated windows use
the descriptive one-sided threshold alpha = 0.01 without a correction across
overlapping windows.

The renderer requires exactly 507 positions-1+2 windows and 252 position-3
windows for each of three branches and two formulas: 4,554 rows in total.
`results/hcv_window_detail.tsv` passed this complete-grid validation. Whole
polyprotein, named-region and codon-position-region results are retained
alongside the compact window tables and publication plots.

## Reproduction

From the repository root:

```bash
python3 research/satute/experiments/005_biological_windows/prepare_hcv_example.py \
  research/satute/artifacts/releases/biological/hcv-codon-windows-20260717/source/hcv_references.gb \
  --outdir WORK/hcv-prepared --write-windows --window-size 36 --window-step 12
```

Infer the full-coding topology and evaluate each window with
`analyze_dataset.py`, using `--model GTR+F+G4` and the frozen topology in
`tree/hcv_polyprotein_codon.treefile`. Then validate and render:

```bash
python3 research/satute/experiments/005_biological_windows/render_hcv_example.py \
  --pos12-detail WORK/pos12/dataset_tree_satute_detail.tsv \
  --pos3-detail WORK/pos3/dataset_tree_satute_detail.tsv \
  --window-manifest WORK/hcv-prepared/window_manifest.tsv \
  --outdir WORK/rendered
```

The frozen release can be re-rendered without repeating the analysis:

```bash
python3 research/satute/experiments/005_biological_windows/render_hcv_example.py \
  --validated-detail research/satute/artifacts/releases/biological/hcv-codon-windows-20260717/results/hcv_window_detail.tsv \
  --outdir research/satute/artifacts/releases/biological/hcv-codon-windows-20260717/results
```

The source GenBank file has SHA-256
`4e9256522f5b86fccd462794204f8f4a6dd54d27e7daba0f2042ff5a053ea0b4`.
All retained payload hashes are listed in `SHA256SUMS`.
