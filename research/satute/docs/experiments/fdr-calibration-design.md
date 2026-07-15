# Branch-level FDR calibration design

The head-to-head benchmark records native Benjamini--Yekutieli decisions for
every pooled branch, but that output is not by itself an empirical FDR study.
The truth of each tested branch must be known independently of the observed
p-values. This document is the gate for any FDR-control claim in the paper.

## Required designs

1. **Exact-null calibration:** 5,000 replicates per tree/model/site cell with
   the focal subtree pair generated under the independence null. Split the
   replicates into 2,500 threshold-training and 2,500 held-out replicates.
2. **Mixed-null operating characteristics:** 2,000 replicates per cell with a
   prespecified mixture of null and finite-branch alternatives. The mixture
   proportions and branch truth labels are fixed before running IQ-TREE.

The same alignment is evaluated by the dominant and eigenvalue-weighted native
implementations. Each replicate must retain a truth-labelled branch map keyed
by normalized split, not by transient branch ID. The audit must contain, for
each formula and decision rule, the number of null branches, alternative
branches, discoveries, false discoveries, true discoveries, FDP and power.

## Acceptance gates

- No truth label may be inferred from the SatuTe p-value or decision.
- Every pooled branch in the FDR family must have exactly one truth label;
  unresolved or duplicated splits fail the run.
- BY decisions must be computed within the prespecified formula-specific tree
  family, with the family size recorded in the audit.
- Report mean FDP (empirical FDR), the distribution of replicate FDPs, power and
  held-out null rejection separately for each formula and correction rule.
- The 5,000/2,000 calibration outputs are separate from the 1,000-replicate
  publication power curves and must not be mixed into their summary tables.

Until these outputs exist and pass the gates, plots labelled `by_fdr` are native
BY decision plots only; they must not be described as validated FDR control.
