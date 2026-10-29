#!/usr/bin/env python3
"""Recompute branch-specific sliding-window SatuTe scores with the current binary.

The curated releases under artifacts/releases/biological/ were produced by a
driver that does not live in this repository, so their window scores could only
be redrawn, never recomputed. This driver closes that gap: it runs SatuTe with
the per-site export, then forms the same linearly weighted rolling z-score the
published analysis used, for every requested branch and both formulas.

Window score, for a window of w sites ending the rolling view:

    weights  v_j = 1 - 0.5 * j / (w - 1)          j = 0 .. w-1
    score    z   = (sum_j v_j c_j / sum_j v_j) / sqrt(satVar / w)

and the score is written at the middle of its window, so a run over n sites
yields n rows of which the first (w-1) - w//2 and the last w//2 are blank.
satVar is read from the pooled .sat.stat row rather than recombined from the
per-category rows, so the mixture variance is the one the test itself used.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

from analyze_dataset import IQTREE_DEFAULT, run_command


def window_weights(window_size):
    if window_size < 2:
        raise ValueError("window size must be at least 2")
    step = 0.5 / (window_size - 1)
    return [1.0 - step * index for index in range(window_size)]


def window_scores(coherence, variance, window_size):
    """Rolling weighted z-scores, centred, with None where the window is short."""
    if variance <= 0.0:
        raise ValueError(f"non-positive satVar: {variance}")
    weights = window_weights(window_size)
    weight_sum = sum(weights)
    scale = math.sqrt(variance / window_size)
    # pandas put the rolling value at the window's last site and then shifted it
    # back by w//2; both steps together place it at first + (w-1) - w//2.
    offset = (window_size - 1) - window_size // 2
    scores = [None] * len(coherence)
    for first in range(len(coherence) - window_size + 1):
        total = 0.0
        for index, weight in enumerate(weights):
            total += weight * coherence[first + index]
        scores[first + offset] = (total / weight_sum) / scale
    return scores


def read_stat(path):
    """Pooled satVar and split membership per (branch id, formula)."""
    branches = {}
    with path.open("r", encoding="utf-8") as handle:
        rows = csv.DictReader(
            (line for line in handle if not line.startswith("#")), delimiter="\t"
        )
        for row in rows:
            if row.get("RateCategory") != "pooled":
                continue
            branch_id = row["ID"]
            entry = branches.setdefault(
                branch_id,
                {
                    "id": branch_id,
                    "label": row.get("Label", ""),
                    "split": tuple(
                        name for name in (row.get("Split") or "").split(",") if name
                    ),
                    "variance": {},
                },
            )
            entry["variance"][row["Formula"]] = float(row["satVar"])
    if not branches:
        raise SystemExit(f"No pooled rows in {path}; was SatuTe run?")
    return branches


def read_sites(path, wanted):
    """Per-site coherence for the requested branch ids, in alignment order."""
    series = {}
    with path.open("r", encoding="utf-8") as handle:
        rows = csv.DictReader(
            (line for line in handle if not line.startswith("#")), delimiter="\t"
        )
        for row in rows:
            key = (row["ID"], row["Formula"])
            if row["ID"] not in wanted:
                continue
            series.setdefault(key, []).append(
                (int(row["Site"]), float(row["SiteCoherence"]))
            )
    if not series:
        raise SystemExit(
            f"{path} holds no rows for branches {sorted(wanted)}; "
            "run with --satute-sites and check the branch selectors"
        )
    ordered = {}
    for key, pairs in series.items():
        pairs.sort()
        ordered[key] = [value for _site, value in pairs]
    return ordered


def resolve_branch(selector, branches):
    """A numeric branch ID, or a taxon name naming that taxon's pendant branch."""
    if selector in branches:
        return selector
    matches = [
        entry["id"]
        for entry in branches.values()
        if entry["split"] == (selector,) or entry["label"] == selector
    ]
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise SystemExit(
            f"No branch matches {selector!r}. Use --list-branches to see the IDs."
        )
    raise SystemExit(f"{selector!r} matches branches {matches}; use a numeric ID.")


def list_branches(branches):
    print(f"{'ID':>6}  {'taxa':>5}  label / split")
    for entry in sorted(branches.values(), key=lambda item: int(item["id"])):
        split = ",".join(entry["split"])
        if len(split) > 90:
            split = split[:87] + "..."
        print(f"{entry['id']:>6}  {len(entry['split']):>5}  {entry['label'] or split}")


def write_window_csv(path, name, scores):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write(f"{name}\n")
        for score in scores:
            handle.write('""\n' if score is None else f"{score!r}\n")


def build_command(args, prefix):
    cmd = [
        str(args.iqtree),
        "-s",
        str(args.alignment),
        "-m",
        args.model,
        "--satute",
        "--satute-sites",
        "--satute-alpha",
        str(args.alpha),
        "--prefix",
        str(prefix),
        "-T",
        str(args.threads),
        "--redo",
    ]
    if args.tree:
        cmd.extend(["-te", str(args.tree)])
        if args.fix_branch_lengths:
            cmd.append("-blfix")
    if args.seed is not None:
        cmd.extend(["--seed", str(args.seed)])
    for extra in args.iqtree_arg:
        cmd.extend(extra.split())
    return cmd


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Recompute SatuTe sliding-window scores for chosen branches."
    )
    parser.add_argument("--alignment", type=Path)
    parser.add_argument("--tree", type=Path, default=None)
    parser.add_argument("--model", default="MFP")
    parser.add_argument(
        "--branch",
        action="append",
        default=[],
        metavar="NAME=SELECTOR",
        help="Output name and branch, by numeric ID or by a taxon on a pendant branch.",
    )
    parser.add_argument("--window-size", type=int, default=36)
    parser.add_argument("--alpha", type=float, default=0.01)
    parser.add_argument("--iqtree", default=str(IQTREE_DEFAULT))
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--dataset", default=None, help="Sub-directory under each formula.")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--fix-branch-lengths", action="store_true")
    parser.add_argument("--iqtree-arg", action="append", default=[])
    parser.add_argument(
        "--reuse-prefix",
        type=Path,
        default=None,
        help="Read an existing .sat.stat/.sat.site pair instead of running IQ-TREE.",
    )
    parser.add_argument("--list-branches", action="store_true")
    args = parser.parse_args(argv)
    if args.reuse_prefix is None and args.alignment is None:
        parser.error("--alignment is required unless --reuse-prefix is given")
    if not args.list_branches and not args.branch:
        parser.error("give at least one --branch NAME=SELECTOR, or --list-branches")
    return args


def main(argv=None):
    args = parse_args(argv)
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    if args.reuse_prefix is not None:
        prefix = args.reuse_prefix
    else:
        prefix = outdir / "satute" / Path(args.alignment).stem
        prefix.parent.mkdir(parents=True, exist_ok=True)
        cmd = build_command(args, prefix)
        run_command(cmd, prefix.parent / "iqtree_command.log")

    branches = read_stat(Path(f"{prefix}.sat.stat"))
    if args.list_branches:
        list_branches(branches)
        return 0

    targets = {}
    for item in args.branch:
        if "=" not in item:
            raise SystemExit(f"--branch wants NAME=SELECTOR, got {item!r}")
        name, selector = item.split("=", 1)
        targets[name] = resolve_branch(selector, branches)

    series = read_sites(Path(f"{prefix}.sat.site"), set(targets.values()))
    dataset = args.dataset or Path(args.alignment or prefix).stem

    written = []
    for name, branch_id in targets.items():
        for formula, variance in sorted(branches[branch_id]["variance"].items()):
            coherence = series.get((branch_id, formula))
            if coherence is None:
                raise SystemExit(f"No per-site rows for branch {branch_id} / {formula}")
            if len(coherence) < args.window_size:
                raise SystemExit(
                    f"{len(coherence)} sites is shorter than the {args.window_size}-site window"
                )
            scores = window_scores(coherence, variance, args.window_size)
            path = (
                outdir
                / formula
                / dataset
                / f"{name}_sliding_window_size_{args.window_size}.csv"
            )
            write_window_csv(path, name, scores)
            finite = sum(1 for score in scores if score is not None)
            written.append((path, finite, len(scores)))

    for path, finite, total in written:
        print(f"{path}  {finite} windows over {total} sites")
    return 0


if __name__ == "__main__":
    sys.exit(main())
