"""Explicit-input workflow adapter for biological sliding windows."""

from satute_analysis.experiments import ExperimentError


def run(context, *, profile, figures):
    del profile, figures
    raise ExperimentError(
        "005_biological_windows requires explicit alignment inputs; "
        f"see {context.experiment.directory / 'README.md'}"
    )


__all__ = ["run"]
