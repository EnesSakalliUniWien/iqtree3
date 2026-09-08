"""Stable project paths independent of the caller's working directory."""

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
REPOSITORY_ROOT = PROJECT_ROOT.parents[1]
