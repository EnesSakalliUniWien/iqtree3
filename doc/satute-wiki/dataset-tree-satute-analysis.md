# Dataset Tree/SatuTe Analysis Driver

`dataset_tree_satute_analysis.py` runs a compact dataset-level workflow:

1. infer or optimize an IQ-TREE maximum-likelihood tree for each alignment,
2. run native `--satute` on the optimized tree,
3. optionally compute neighborhood-style branch support labels with SH-aLRT,
   UFBoot, or both,
4. join tree likelihood, support labels, and SatuTe saturation rows by branch
   split.

Example:

```bash
python3 doc/satute-wiki/dataset_tree_satute_analysis.py \
  example/example.phy \
  --iqtree ./build/iqtree3 \
  --outdir /tmp/iqtree-satute-dataset-analysis \
  --model JC \
  --support-mode alrt \
  --support-replicates 1000 \
  --threads 1 \
  --quiet
```

For exploratory runs without branch-support resampling:

```bash
python3 doc/satute-wiki/dataset_tree_satute_analysis.py \
  /path/to/alignment.fa \
  --iqtree ./build/iqtree3 \
  --outdir /tmp/iqtree-satute-dataset-analysis \
  --model MFP \
  --support-mode none
```

Outputs:

- `dataset_tree_satute_detail.tsv`: one pooled SatuTe row per branch and
  formula, plus tree log-likelihood, AIC/BIC, branch split, SH-aLRT/UFBoot
  support fields when available, and saturation decisions.
- `dataset_tree_satute_summary.tsv`: per-dataset/per-formula counts of
  informative versus saturated branches, Bonferroni decisions, median `satZ`
  and `satP`, high-support saturation counts, and likelihood scores.
- `manifest.json`: exact IQ-TREE commands used for the run.
- `runs/<dataset>/`: native IQ-TREE outputs, including `.iqtree`, `.treefile`,
  `.sat.stat`, `.sat.tree`, and command log.

The support join is split-based. A row can have `tree_split_joined=1` while
`support_present=0` for terminal branches or for runs without support labels.
For `--support-mode both`, IQ-TREE usually writes internal labels as
`SH-aLRT/UFBoot`; for `alrt` or `ufboot` alone, the single numeric label is
reported in the matching `sh_alrt` or `ufboot` column.
