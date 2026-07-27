# 005 - Biological sliding windows

Contains the curated and independent biological window analyses used as
examples in the manuscript. These are distinct from the paired simulation
power experiments.

- Dataset analysis: `analyze_dataset.py`
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
