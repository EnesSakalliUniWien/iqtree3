# Paper Map

The paper introduces SatuTe as a branch-level saturation test. The central
question is not whether a whole alignment has signal in the abstract, but
whether the two subtrees separated by a particular branch still share enough
phylogenetic information to justify that branch. This shifts saturation from a
pairwise sequence problem to a subtree-to-subtree problem.

The paper was published in Molecular Biology and Evolution as volume 42, issue
5, article `msaf090`, with DOI `10.1093/molbev/msaf090`. The official article
page lists the publication date as 27 May 2025 and links the implementation and
auxiliary analysis scripts in the data availability section.

The paper states that the implementation is available at:

```text
https://github.com/Elli-ellgard/SatuTe
```

and that auxiliary scripts are available at:

```text
https://github.com/Elli-ellgard/SatuTe-example-analyses
```

The example-analysis repository contains the downstream Tree-of-Life analysis
workflow and data used for the biological examples. The paper's main simulation
study is described in the article and supplementary methods rather than as a
single executable reproduction script in the example repository. The notes
therefore separate the published simulation target from the smaller development
script used to exercise the IQ-TREE integration.

## Fig. 2 Simulation Design

The simulation section evaluates SatuTe on two trees. The first is a five-taxon
tree where the branch of interest is external. The second is a balanced
16-taxon tree where the branch of interest is the central internal branch
separating two eight-taxon subtrees.

The paper varies the focal branch length over:

```text
0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1.0, 1.5,
2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.5, 10.0
```

It uses alignments with:

```text
100, 1000, and 10000 sites
```

and simulates 1000 DNA alignments for each parameter combination under the
Jukes-Cantor model. Each alignment is analyzed in four scenarios: the true tree
with fixed branch lengths, the true topology with ML-estimated branch lengths,
the ML-inferred tree, and the ML-inferred tree with Bonferroni correction. The
reported curve is the fraction of replicates in which the target branch is
called phylogenetically informative.

## Development Smoke Test

The script [reconstruct_fig2_simulations.py](reconstruct_fig2_simulations.py)
implements the same four scenario labels using IQ-TREE AliSim and the native
`--satute` integration. This is useful for testing branch selection and output
parsing, but it is not yet a paper-faithful reproduction because the original
study used Seq-Gen and EvoNAPS-derived 16-taxon subtree branch lengths.

The first reduced run used two replicates, site lengths of 100 and 1000, and
branch lengths of 0.1, 1.0, and 5.0. This design is much smaller than the
published simulation grid and uses different simulation details. It is a
development check for command-line wiring, not an estimate of the published
power curves.
