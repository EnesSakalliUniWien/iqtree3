# SatuTe manuscript editorial log

## Deliverables

- `UPDATED_main.tex`: clean LaTeX entry point.
- `UPDATED_sections/`: edited narrative and methods sections.
- `UPDATED_main.pdf`: compiled 18-page manuscript.

The existing `main.tex`, original section files, and `manuscript.pdf` were left unchanged.

## Structural revisions

1. Reframed original SatuTe as the primary model in the title, abstract, introduction, Results, Methods, and Discussion.
2. Moved the original JC benchmark ahead of the model-extension benchmarks.
3. Presented eigenvalue weighting as a model-aware extension and sensitivity analysis rather than as a replacement for original SatuTe.
4. Separated the invariant-site correction from the primary SatuTe formulation.
5. Divided the empirical Results into Tree-of-Life and HCV applications, each led by the original SatuTe finding and followed by the weighted sensitivity result.
6. Added explicit limits on interpretation: SatuTe localizes signal but does not select a true topology, establish model adequacy, or infer viral population structure.

## Language and style revisions

- Tightened the abstract around the problem, primary test, extensions, principal results, and conclusion.
- Removed repetition between the Introduction and Discussion.
- Shortened long sentences and clarified antecedents.
- Standardized the clean manuscript to American English.
- Replaced ambiguous uses of “published” and “dominant” with “original SatuTe” when referring to the primary statistic.
- Preserved all reported numerical results, equations, citations, figure files, and substantive claims.

## Validation

- Academic Writing Toolkit paragraph-logic review: no issues in the original review excerpts and no issues in the revised abstract and discussion.
- LaTeX compilation: successful with TeX Live/`latexmk`.
- Cross-references and citations: no undefined references or citations.
- Visual check: title page, equations, captions, and full-page workflow figure render correctly.
