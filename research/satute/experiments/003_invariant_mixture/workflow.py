"""Smoke, full-run, and rendering adapter for invariant mixtures."""


def run(context, *, profile, figures):
    output = context.output_dir(profile)
    if figures:
        context.run_python(
            "render_full.py",
            "--input-dir",
            output,
            "--output-prefix",
            output / "figures/full-analysis",
            "--png-scale",
            3,
        )
        return

    arguments = [
        "--outdir",
        output,
        "--figure-prefix",
        output / "figures/overview",
        "--png-scale",
        3,
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
        arguments.extend(("--null-reps", 20000, "--alternative-reps", 20000))
    context.run_python("run.py", *arguments)


__all__ = ["run"]
