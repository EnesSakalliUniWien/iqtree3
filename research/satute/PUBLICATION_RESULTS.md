# SatuTe publication-suite results

This index points to the verified local outputs generated from the completed
LiSC publication and calibration arrays. The manuscript-facing figures use
fixed GTR+F, GTR+F+G4, GTR+F+I+G4 and fixed protein +G4 models. Each
publication cell contains 1,000 paired replicates.

## Built paper and overview curves

- [Built manuscript](manuscript/manuscript.pdf)
- [Exploratory power-curve summary](manuscript/figures/figure_power_curve_summary.pdf)
- [Dense invariant-mixture power analysis](manuscript/figures/figure_weighted_invariant_full.pdf)
- [Tracked publication figure subset](manuscript/figures/publication-suite/README.md)

## Matched nucleotide suite

The largest weighted-minus-dominant informative-call increase for matched
GTR+F+G4 was 0.241 (95% interval 0.215--0.269) without adjustment, 0.576
(0.546--0.607) under taxon Bonferroni, and 0.414 (0.384--0.446) under
Benjamini--Yekutieli. The full figure tree contains matched GTR+F, GTR+F+G4
and GTR+F+I+G4 model pairs under all three rules:

- [Matched GTR figure tree](artifacts/cluster/lisc/paper-matched-gtr-rate-invar-reps1000-ff7c9034b00-20260715/figures-r)
- [Matched GTR+G4 unadjusted paired curve](artifacts/cluster/lisc/paper-matched-gtr-rate-invar-reps1000-ff7c9034b00-20260715/figures-r/unadjusted/GTR_PF06346_G4__GTR_PF06346_G4/figure_paired_difference_2d_rule_unadjusted__sim_GTR_PF06346_G4__eval_GTR_PF06346_G4.pdf)
- [Matched GTR+I+G4 taxon-Bonferroni paired curve](artifacts/cluster/lisc/paper-matched-gtr-rate-invar-reps1000-ff7c9034b00-20260715/figures-r/taxon_bonferroni/GTR_PF06346_I_G4__GTR_PF06346_I_G4/figure_paired_difference_2d_rule_taxon_bonferroni__sim_GTR_PF06346_I_G4__eval_GTR_PF06346_I_G4.pdf)

## Protein suite

All 1,728 protein model--rule--grid cells were nonnegative: 831 had a positive
weighted-minus-dominant difference and 897 were equal. The largest differences
were as follows.

| Model | Unadjusted | Taxon Bonferroni | Benjamini--Yekutieli |
| --- | ---: | ---: | ---: |
| LG+G4 | 0.030 (0.020--0.041) | 0.327 (0.298--0.356) | 0.094 (0.076--0.113) |
| WAG+G4 | 0.075 (0.059--0.092) | 0.469 (0.438--0.499) | 0.165 (0.142--0.188) |
| JTT+G4 | 0.479 (0.447--0.510) | 0.555 (0.524--0.587) | 0.514 (0.483--0.545) |
| Q.pfam+G4 | 0.041 (0.029--0.054) | 0.354 (0.325--0.384) | 0.118 (0.099--0.138) |

- [Complete protein figure tree](artifacts/cluster/lisc/paper-matched-protein-g4-reps1000-ff7c9034b00-20260715/figures-r)
- [JTT+G4 taxon-Bonferroni paired curve](artifacts/cluster/lisc/paper-matched-protein-g4-reps1000-ff7c9034b00-20260715/figures-r/taxon_bonferroni/JTT_G4__JTT_G4/figure_paired_difference_2d_rule_taxon_bonferroni__sim_JTT_G4__eval_JTT_G4.pdf)

## GTR model misspecification

For GTR+F simulations evaluated under JC, K2P and F81, the largest paired
increase was 0.131 (0.107--0.157) without adjustment, 0.206
(0.176--0.236) under taxon Bonferroni, and 0.209 (0.180--0.240) under
Benjamini--Yekutieli. The maxima occurred in the K2P evaluation.

- [Complete misspecification figure tree](artifacts/cluster/lisc/paper-gtr-misspec-reps1000-ff7c9034b00-20260715/figures-r)
- [K2P Benjamini--Yekutieli paired curve](artifacts/cluster/lisc/paper-gtr-misspec-reps1000-ff7c9034b00-20260715/figures-r/by_fdr/GTR_PF06346__K2P/figure_paired_difference_2d_rule_by_fdr__sim_GTR_PF06346__eval_K2P.pdf)

## Truth-labelled FDR and power calibration

The calibration design used 5,000 exact-null replicates per cell, split into
2,500 training and 2,500 held-out replicates, plus an independent 2,000-replicate
50:50 mixed-null design. In the GTR mixed-null grid, power ranged from 0.929 to
0.989 and empirical FDR from zero to 0.00386. Formula-paired power differences
were between -0.000071 and 0.000571. Corrected tie-aware training thresholds
gave held-out null rejection rates from 0.0108 to 0.0644 across all 96 cells;
no cell retained the erroneous 100% rejection produced by splitting a discrete
adjusted-p-value tie block.

- [Corrected calibration figure tree](artifacts/cluster/lisc/fdr-calibration-reps5000-2000-corrected-ties/figures)
- [Weighted GTR FDR and power](artifacts/cluster/lisc/fdr-calibration-reps5000-2000-corrected-ties/figures/figure_eigenvalue_weighted_fdr_power_2d_model_GTR_PF06346.pdf)
- [Weighted GTR held-out calibration](artifacts/cluster/lisc/fdr-calibration-reps5000-2000-corrected-ties/figures/figure_eigenvalue_weighted_null_calibration_2d_model_GTR_PF06346.pdf)
- [Operating summary table](artifacts/cluster/lisc/fdr-calibration-reps5000-2000-corrected-ties/figures/fdr_operating_summary.tsv)
- [Paired-difference summary table](artifacts/cluster/lisc/fdr-calibration-reps5000-2000-corrected-ties/figures/fdr_paired_difference_summary.tsv)
- [Corrected held-out calibration table](artifacts/cluster/lisc/fdr-exact-null-reps5000-a932abf97b2-20260715-postprocessed-corrected-ties/heldout_null_calibration.tsv)

## Figure QA

The matched nucleotide, matched protein and misspecification trees contain,
respectively, 144, 192 and 144 nonempty outputs. Every model pair and decision
rule has independent dominant, eigenvalue-weighted and paired-difference 2D
PDF/SVG/450-dpi PNG outputs, plus vector 3D PDF and self-contained HTML. The
corrected calibration bundle contains 22 PDFs, 10 SVGs, 10 450-dpi PNGs, 12
self-contained HTML files and three summary tables. All PDFs and SVGs parse,
and no external HTML dependency directory remains.
