"""Smoke, full-run, and rendering adapter for relative weighting."""


def run(context, *, profile, figures):
    output = context.output_dir(profile)
    if figures:
        arguments = [
            "--summary",
            output / "head_to_head_summary.tsv",
            "--outdir",
            output / "figures",
            "--expected-reps",
            1 if profile == "smoke" else 1000,
            "--tree-cases",
            "five_external,sixteen_internal",
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
        "--simulation-models",
        "JC",
        "--evaluation-models",
        "JC",
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
                "--evonaps-branch-lengths",
                context.project_root
                / "references/evonaps/evonaps_16taxon_branch_lengths.tsv",
                "--resume",
            )
        )
    else:
        arguments.extend(
            (
                "--reps",
                1000,
                "--paper-grid",
                "--tree-cases",
                "five_external,sixteen_internal",
                "--evonaps-branch-lengths",
                context.project_root
                / "references/evonaps/evonaps_16taxon_branch_lengths.tsv",
            )
        )
    context.run_python("run.py", *arguments)


__all__ = ["run"]
