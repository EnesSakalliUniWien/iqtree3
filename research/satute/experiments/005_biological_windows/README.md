# 005 - Biological sliding windows

Contains the curated and independent biological window analyses used as
examples in the manuscript. These are distinct from the paired simulation
power experiments.

- Dataset analysis: `analyze_dataset.py`
- Sliding-window recomputation: `sliding_window.py`
- Curated renderer: `render_curated.py`
- Independent renderer: `render_independent.py`
- HCV reference preparation: `prepare_hcv_example.py`
- HCV grid validator and renderer: `render_hcv_example.py`
- Manifest and explicit-input instructions: `experiment.json`, this README
- Reviewed inputs/results: `../../artifacts/releases/biological/`

The HCV example uses ten accession-versioned genotype references, a
codon-preserving polyprotein alignment and equal 36-site windows from codon
positions 1+2 or position 3. Its complete 4,554-row target grid, figures and
provenance are frozen in
`../../artifacts/releases/biological/hcv-codon-windows-20260717/`.

## Recomputing the sliding windows

`sliding_window.py` recomputes the branch-specific window scores that the
curated releases only stored. It runs SatuTe with `--satute-sites`, reads the
per-site coherence and the pooled `satVar`, and writes one CSV per branch and
formula in the release layout `<formula>/<dataset>/<name>_sliding_window_size_<w>.csv`.

The published analyses used 36-site windows at alpha 0.01, under `LG+G4` for
`protein_based_2D_tree` and `GTR+G4` for `rRNA_based_3D_tree`, on the branches
leading to Eukaryota and to yeast.

Both Tree-of-Life alignments carry thousands of taxa, so ask for the branch
table first and then restrict the scoring run to the branches you need:

```bash
experiments/005_biological_windows/sliding_window.py \
  --alignment <data>/rRNA_based_3D_tree/rRNA_alignment.fasta \
  --tree <data>/rRNA_based_3D_tree/rRNA_based_3D_ToL.treefile \
  --fix-branch-lengths --model GTR+G4 --alpha 0.01 \
  --outdir artifacts/local/biological-windows/rRNA --list-branches
```

```bash
experiments/005_biological_windows/sliding_window.py \
  --alignment <data>/rRNA_based_3D_tree/rRNA_alignment.fasta \
  --tree <data>/rRNA_based_3D_tree/rRNA_based_3D_ToL.treefile \
  --fix-branch-lengths --model GTR+G4 --alpha 0.01 \
  --iqtree-arg "--satute-edges edges.txt" \
  --dataset rRNA_based_3D_tree \
  --branch branch_to_eukaryota=<id> --branch branch_to_yeast=<taxon> \
  --outdir artifacts/local/biological-windows/rRNA
```

A branch selector is either a numeric ID from `--list-branches` or, for a
pendant branch, the taxon name itself. `--reuse-prefix` skips the IQ-TREE run
and reads an existing `.sat.stat`/`.sat.site` pair.

The scores are not comparable to the frozen June 2026 releases: those predate
the centred mixture statistic, so a rerun measures the current test, not the
published one.
