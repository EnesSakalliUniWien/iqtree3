"""Stable project paths independent of the caller's working directory."""

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
REPOSITORY_ROOT = PROJECT_ROOT.parents[1]
ARTIFACT_ROOT = PROJECT_ROOT / "artifacts"
LOCAL_ARTIFACT_ROOT = ARTIFACT_ROOT / "local"
RELEASE_ROOT = ARTIFACT_ROOT / "releases"
MANUSCRIPT_ROOT = PROJECT_ROOT / "manuscript"


def relative_to_project(path):
    """Resolve a path against the SatuTe research root."""
    value = Path(path)
    return value if value.is_absolute() else PROJECT_ROOT / value
