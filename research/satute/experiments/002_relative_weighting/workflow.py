"""Smoke, full-run, and rendering adapter for relative weighting."""

from satute_analysis.experiments import ExperimentError


def _csv(values):
    return ",".join(str(value) for value in values)


def _run_fdr_calibration(context, profile, figures, study, configuration):
    if figures:
        raise ExperimentError(
            "FDR calibration figures require merged exact-null, mixed-null, and "
            "held-out tables; use render_fdr_calibration.R after postprocessing"
        )

    design = configuration["fdr_calibration_design"]
    exact = study == "fdr_exact_null"
    if profile == "smoke":
        arguments = [
            "--design",
            "exact_null" if exact else "mixed_null",
            "--reps",
            2,
            "--training-reps",
            1 if exact else 0,
            "--models",
            design["models"][0],
            "--site-lengths",
            design["site_lengths"][0],
            "--branch-lengths",
            design["branch_lengths"][0],
            "--tree-cases",
            design["tree_cases"][0],
        ]
    else:
        arguments = [
            "--design",
            "exact_null" if exact else "mixed_null",
            "--reps",
            (
                design["exact_null_replicates_per_cell"]
                if exact
                else design["mixed_null_replicates_per_cell"]
            ),
            "--training-reps",
            design["exact_null_training_replicates"] if exact else 0,
            "--models",
            _csv(design["models"]),
            "--site-lengths",
            _csv(design["site_lengths"]),
            "--branch-lengths",
            _csv(design["branch_lengths"]),
            "--tree-cases",
            _csv(design["tree_cases"]),
        ]
    context.run_python(
        "fdr_calibration_run.py",
        "--iqtree",
        context.iqtree,
        "--outdir",
        context.output_dir(profile),
        "--evonaps-branch-lengths",
        context.project_root
        / "references/evonaps/evonaps_16taxon_branch_lengths.tsv",
        *arguments,
    )


def run(context, *, profile, figures, study):
    output = context.output_dir(profile)
    configuration = context.experiment.configuration
    if study in {"fdr_exact_null", "fdr_mixed_null"}:
        _run_fdr_calibration(context, profile, figures, study, configuration)
        return

    design = configuration["publication_curve_design"]
    study_design = design["studies"][study]
    if figures:
        arguments = [
            "--summary",
            output / "head_to_head_summary.tsv",
            "--outdir",
            output / "figures",
            "--expected-reps",
            1 if profile == "smoke" else design["replicates_per_cell"],
            "--tree-cases",
            _csv(design["tree_cases"]),
        ]
        if profile == "smoke":
            arguments.append("--allow-incomplete")
        context.run_python("render.py", *arguments)
        return

    arguments = [
        "--iqtree",
        context.iqtree,
        "--outdir",
        output,
    ]
    if profile == "smoke":
        arguments.extend(
            (
                "--reps",
                1,
                "--site-lengths",
                "100",
                "--branch-lengths",
                "0.5,1.0",
                "--tree-cases",
                design["tree_cases"][0],
                "--model-pairs",
                study_design["model_pairs"][0],
                "--scenario-set",
                study_design["scenario_set"],
                "--simulator",
                design["simulator"],
                "--evonaps-branch-lengths",
                context.project_root
                / "references/evonaps/evonaps_16taxon_branch_lengths.tsv",
            )
        )
    else:
        arguments.extend(
            (
                "--reps",
                design["replicates_per_cell"],
                "--site-lengths",
                _csv(design["site_lengths"]),
                "--branch-lengths",
                _csv(design["branch_lengths"]),
                "--tree-cases",
                _csv(design["tree_cases"]),
                "--model-pairs",
                _csv(study_design["model_pairs"]),
                "--scenario-set",
                study_design["scenario_set"],
                "--simulator",
                design["simulator"],
                "--evonaps-branch-lengths",
                context.project_root
                / "references/evonaps/evonaps_16taxon_branch_lengths.tsv",
            )
        )
    context.run_python("run.py", *arguments)


__all__ = ["run"]
