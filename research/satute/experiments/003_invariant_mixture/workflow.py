"""Smoke, full-run, and rendering adapter for invariant mixtures."""


def _csv(values):
    return ",".join(str(value) for value in values)


def run(context, *, profile, figures, study):
    del study
    output = context.output_dir(profile)
    configuration = context.experiment.configuration
    if figures:
        context.run_python(
            "render_full.py",
            "--input-dir",
            output,
            "--output-prefix",
            output / "figures/full-analysis",
            "--png-scale",
            configuration["png_scale"],
        )
        return

    arguments = [
        "--outdir",
        output,
        "--figure-prefix",
        output / "figures/overview",
        "--png-scale",
        configuration["png_scale"],
    ]
    if profile == "smoke":
        arguments.extend(
            (
                "--site-lengths",
                "100,1000",
                "--branch-lengths",
                "20,50,10000",
                "--null-reps",
                200,
                "--alternative-reps",
                200,
            )
        )
    else:
        arguments.extend(
            (
                "--site-lengths",
                _csv(configuration["site_lengths"]),
                "--branch-lengths",
                _csv(configuration["branch_lengths"]),
                "--null-reps",
                configuration["null_replicates"],
                "--alternative-reps",
                configuration["alternative_replicates"],
                "--seed",
                configuration["seed"],
            )
        )
    context.run_python("run.py", *arguments)


__all__ = ["run"]
