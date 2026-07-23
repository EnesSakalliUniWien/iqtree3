#!/usr/bin/env python3
"""Run a registered SatuTe experiment through its local workflow adapter."""

import argparse
import sys
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT / "src"))

from satute_analysis.experiments import ExperimentCatalog, ExperimentError


def main() -> int:
    catalog = ExperimentCatalog.discover(PROJECT)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiment", choices=catalog.ids, required=True)
    parser.add_argument("--profile", choices=("smoke", "full"), default="smoke")
    parser.add_argument("--figures", action="store_true")
    args = parser.parse_args()
    try:
        catalog.run(
            args.experiment,
            profile=args.profile,
            figures=args.figures,
            python=sys.executable,
        )
    except ExperimentError as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    sys.exit(main())
