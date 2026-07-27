"""Smoke, full-run, and rendering adapter for paper reconstruction."""


def run(context, *, profile, figures, study):
    del study
    output = context.output_dir(profile)
    if figures:
        context.run_python(
            "render.py",
            "--summary",
            output / "fig2_reconstruction_summary.tsv",
            "--output",
            output / "figures/fig2_reconstruction.svg",
        )
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
                "five_external,sixteen_internal",
            )
        )
    else:
        arguments.extend(("--reps", 1000, "--paper-grid"))
    context.run_python("run.py", *arguments)


__all__ = ["run"]
