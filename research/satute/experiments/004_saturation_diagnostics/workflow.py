"""Smoke, full-run, and rendering adapter for saturation diagnostics."""


def run(context, *, profile, figures, study):
    del study
    output = context.output_dir(profile)
    configuration = context.experiment.configuration
    arguments = [
        "--iqtree",
        context.iqtree,
        "--outdir",
        output,
        "--figure-prefix",
        output / "figures/main",
        "--category-figure-prefix",
        output / "figures/categories",
    ]
    if profile == "smoke":
        arguments.extend(
            (
                "--sites",
                100,
                "--branch-lengths",
                "0.5,10",
                "--null-length",
                10,
                "--reps",
                2,
                "--null-reps",
                10,
                "--workers",
                2,
            )
        )
    else:
        arguments.extend(
            (
                "--sites",
                configuration["sites"],
                "--reps",
                configuration["alternative_replicates"],
                "--null-reps",
                configuration["null_replicates"],
            )
        )
    if figures:
        arguments.append("--plot-only")
    context.run_python("run.py", *arguments)


__all__ = ["run"]
