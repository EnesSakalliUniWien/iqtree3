# SatuTe — branch saturation test

Tests whether a branch still carries detectable phylogenetic signal between the
two subtrees it separates, or whether substitution has saturated it. Native
IQ-TREE implementation of:

> Manuel C, Sakalli E, Schmidt HA, Viñas C, von Haeseler A, Elgert C.
> When the Past Fades: Detecting Phylogenetic Signal with SatuTe.
> *Molecular Biology and Evolution* 42(5):msaf090, 2025.
> [doi:10.1093/molbev/msaf090](https://doi.org/10.1093/molbev/msaf090)

## Build

Standard IQ-TREE build; no extra dependencies.

```bash
cmake -S . -B build && cmake --build build --target iqtree3 -j
```

## Run

```bash
iqtree3 -s aln.fa -te tree.nwk -m GTR+F -blfix --satute
```

| flag | meaning |
|---|---|
| `--satute` | run the analysis |
| `--satute-alpha NUM` | significance and FDR level (default `0.05`) |
| `--satute-edges FILE` | restrict to branch IDs listed in `FILE` |
| `--satute-sites` | also write per-site contributions (`.sat.site`) |
| `--satute-sites-max NUM` | cap `.sat.site` to NUM branches, highest `satP` first |
| `--satute-enhanced-weights` | experimental GLS mode weighting; off by default |

Writes `<prefix>.sat.stat` (one row per branch × formula × rate category),
plus `.sat.branch`, `.sat.tree` and `.sat.tree.nex`. Key columns: `satC`
(coherence), `satZ`, `satP`, `Decision`, and `DecisionFDR`.

Requires a reversible model, an unrooted single-partition tree, and no
ascertainment-bias correction; mixture models are not supported.

## Two formulas

Every branch is reported twice.

- **`dominant`** — the published statistic. Uses only the eigenmodes tied to
  the largest non-zero eigenvalue, selected by a hard tolerance.
- **`eigenvalue_weighted`** — uses every non-zero mode, weighting mode *k* by
  `exp((λ_k − λ_max)·t)`. No selection step, so it cannot flip its mode set
  between replicates when eigenvalues are near-degenerate.

The statistic is **centred**: it accumulates a covariance
`Σ_k w_k (E[L_k R_k] − E[L_k]E[R_k])`, not the raw cross-moment. The uncentred
form is only a valid independence test when `E[L_k]=E[R_k]=0`, which fails
under a misspecified model and biases the z-score by a factor that grows with
alignment length.

BY-FDR is applied separately within each formula's pooled branch rows.
Rate-category rows are diagnostic and are not members of either FDR family.

## Layout

| file | role |
|---|---|
| `../satute.cpp` | driver: walks branches, builds weights, pools, writes |
| `satute_statistics.cpp` | the kernel — partials, coherence, variance, z |
| `satute_formulas.cpp` | eigenvalue selection and mode weights |
| `satute_rates.cpp` | rate-category extraction |
| `satute_information.cpp` | information fraction / saturation index |
| `satute_topology.cpp` | splits and taxon counts |
| `satute_edges.cpp` | `--satute-edges` parsing |
| `satute_results.cpp` | row construction and pooling |
| `satute_report_writer.cpp` | output files |
| `satute_format.cpp` | number and vector formatting |
| `satute_types.h` | shared structs |

## Performance

The pattern loop is OpenMP-parallel (`-T`), with per-thread accumulators
combined in fixed order so results are bit-reproducible for a given thread
count.

Every branch needs the two directional partial likelihoods that meet on it.
Computing those per branch re-walks a whole subtree for each of the `2n−3`
branches, which is `O(n²·patterns)`. Instead all `2(2n−3)` directional messages
are built once per rate category by a two-pass sweep (`SatuTeDirectionalCache`),
making the statistic `O(n·patterns)`.

Measured on an M-series Mac, GTR, 2000 sites, `-T 1`:

| taxa | before | after | speedup |
|---|---|---|---|
| 20 | 0.08 s | 0.02 s | 4× |
| 50 | 0.49 s | 0.03 s | 16× |
| 100 | 2.09 s | 0.06 s | 35× |
| 200 | 10.91 s | 0.11 s | 99× |

Linear in alignment length throughout (16 000 sites at 50 taxa: 3.92 s → 0.21 s).

The cache holds `nodes × max_degree × patterns × states` doubles. Rate category
is the outer loop, so each category's cache is released before the next is
built and only one is resident: peak memory is one category's worth regardless
of how many the model has. If that single cache would exceed 2 GB, SatuTe says
so and falls back to the per-branch recursion, which is slower but produces
identical numbers.

Setting `SATUTE_DISABLE_DIRECTIONAL_CACHE=1` forces that fallback, so the two
paths can be diffed against each other; they are verified byte-identical across
JC, GTR, +G4, +I+G4, +R4, protein, `--satute-sites` and
`--satute-enhanced-weights`.
