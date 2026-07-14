"""Shared mathematical and project utilities for SatuTe experiments."""

from .project import PROJECT_ROOT, REPOSITORY_ROOT

__all__ = [
    "PROJECT_ROOT",
    "REPOSITORY_ROOT",
    "build_q",
    "parse_model",
    "reversible_eigendecomposition",
    "transition_matrix",
]


def __getattr__(name):
    """Load NumPy-backed model helpers only when they are requested."""

    if name in {"build_q", "parse_model", "reversible_eigendecomposition", "transition_matrix"}:
        from . import models

        return getattr(models, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
