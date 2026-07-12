#!/usr/bin/env python3
"""Plot the curated independent biological sliding-window comparison."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from render_curated import (
    convert_svg,
    line,
    nice_ticks,
    polyline,
    read_window_scores,
    rect,
    svg_document,
    text,
)


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DEFAULT_RUN_DIR = (
    REPO_ROOT
    / "artifacts"
    / "releases"
    / "biological"
    / "independent_enhanced_sliding_window_20260619_151135"
)

DATASET_LABELS = {
    "protein_based_2D_tree": "protein 2D tree",
    "rRNA_based_3D_tree": "16S rRNA 3D tree",
}
BRANCH_LABELS = {
    "branch_to_eukaryota": "Eukaryota",
    "branch_to_yeast": "yeast",
}
SERIES = (
    ("published", "published SatuTe", "#0A3FBA", "published_satute_csv", "published_saturated_windows", "published_saturated_fraction"),
    (
        "eigenvalue_weighted",
        "eigenvalue-weighted",
        "#E69504",
        "eigenvalue_weighted_csv",
        "eigenvalue_weighted_saturated_windows",
        "eigenvalue_weighted_saturated_fraction",
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot curated independent biological SatuTe comparisons."
    )
    parser.add_argument(
        "--run-dir",
        type=Path,
        default=DEFAULT_RUN_DIR,
        help="Curated independent biological rerun directory.",
    )
    return parser.parse_args()


def resolve_path(path_text: str) -> Path:
    path = Path(path_text)
    if path.is_absolute():
        return path
    for candidate in (REPO_ROOT / path, Path.cwd() / path):
        if candidate.exists():
            return candidate
    raise FileNotFoundError(path_text)


def read_summary(run_dir: Path) -> list[dict[str, str]]:
    summary = run_dir / "comparison_to_manuscript" / "published_vs_eigenvalue_weighted_summary.tsv"
    with summary.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    expected = {
        (dataset, branch)
        for dataset in DATASET_LABELS
        for branch in BRANCH_LABELS
    }
    observed = {(row["dataset"], row["branch_id"]) for row in rows}
    missing = sorted(expected - observed)
    extra = sorted(observed - expected)
    if missing or extra:
        raise ValueError(f"Unexpected dataset or branch set: missing={missing}, extra={extra}")
    return rows


def panel_title(row: dict[str, str]) -> str:
    return f"{DATASET_LABELS[row['dataset']]} - {BRANCH_LABELS[row['branch_id']]}"


def plot_curves(rows: list[dict[str, str]], run_dir: Path) -> None:
    width, height = 1550, 1180
    plot_left, plot_right = 100, width - 35
    plot_width = plot_right - plot_left
    panel_height = 205
    panel_tops = [80, 355, 630, 905]
    body: list[str] = []

    legend_y = 34
    legend_x = 1060
    for index, (_, label, color, _, _, _) in enumerate(SERIES):
        x = legend_x + index * 220
        body.append(line(x, legend_y - 6, x + 58, legend_y - 6, color=color, width=3.0))
        body.append(text(x + 72, legend_y, label, size=18, anchor="start"))

    for row, top in zip(rows, panel_tops, strict=True):
        windows = int(row["windows"])
        series_values = [
            read_window_scores(resolve_path(row[csv_key]), windows)
            for _, _, _, csv_key, _, _ in SERIES
        ]
        y_min = min(-1.2, *(min(values) - 0.4 for values in series_values))
        y_max = max(max(values) for values in series_values) + 1.0
        ticks = nice_ticks(y_min, y_max, 5.0 if y_max > 20 else 2.5)
        bottom = top + panel_height

        def sx(index: int) -> float:
            return plot_left + plot_width * index / max(1, windows - 1)

        def sy(value: float) -> float:
            return bottom - panel_height * (value - y_min) / (y_max - y_min)

        body.append(text(width / 2, top - 18, panel_title(row), size=20, weight="700"))
        body.append(rect(plot_left, top, plot_width, panel_height, fill="#ffffff", stroke="#111111"))

        for tick in ticks:
            if y_min <= tick <= y_max:
                y = sy(tick)
                body.append(line(plot_left, y, plot_right, y, color="#d7d7d7", width=0.8))
                body.append(text(plot_left - 14, y + 5, f"{tick:g}", size=14, anchor="end"))
        for tick in range(0, windows, 500):
            x = sx(tick)
            body.append(line(x, top, x, bottom, color="#eeeeee", width=0.8))
            if top == panel_tops[-1]:
                body.append(text(x, bottom + 22, str(tick), size=14))

        body.append(line(plot_left, sy(2.326347874388028), plot_right, sy(2.326347874388028), color="#999999", width=1.3, dash="8 6"))
        for values, (_, _, color, _, _, _) in zip(series_values, SERIES, strict=True):
            body.append(polyline([(sx(i), sy(value)) for i, value in enumerate(values)], color=color, width=1.3))
        body.append(text(34, top + panel_height / 2, "Window z-score", size=16, rotate=-90))

    body.append(text(width / 2, height - 16, "Window index", size=18))
    out_base = run_dir / "comparison_to_manuscript" / "published_vs_eigenvalue_weighted_curves"
    convert_svg(svg_document(width, height, body), out_base, width)


def plot_saturated_windows(rows: list[dict[str, str]], run_dir: Path) -> None:
    width, height = 1120, 760
    panel_width, panel_height = 420, 235
    panel_positions = [(95, 80), (620, 80), (95, 455), (620, 455)]
    body: list[str] = []

    for row, (left, top) in zip(rows, panel_positions, strict=True):
        bottom = top + panel_height
        values = [100.0 * float(row[fraction_key]) for *_, fraction_key in SERIES]
        y_max = max(6.0, max(values) * 1.25)
        step = 10.0 if y_max > 20 else 1.0
        body.append(text(left + panel_width / 2, top - 26, panel_title(row), size=20, weight="700"))
        body.append(rect(left, top, panel_width, panel_height, fill="#ffffff", stroke="#111111"))

        for tick in nice_ticks(0, y_max, step):
            if 0 <= tick <= y_max:
                y = bottom - panel_height * tick / y_max
                body.append(line(left, y, left + panel_width, y, color="#d7d7d7", width=0.8))
                body.append(text(left - 12, y + 5, f"{tick:g}", size=13, anchor="end"))

        bar_width = 115
        xs = [left + 70, left + 235]
        for x, value, (_, label, color, _, count_key, _) in zip(xs, values, SERIES, strict=True):
            bar_height = panel_height * value / y_max
            body.append(rect(x, bottom - bar_height, bar_width, bar_height, fill=color))
            body.append(text(x + bar_width / 2, bottom - bar_height - 8, row[count_key], size=15))
            for index, part in enumerate(label.split(" ")):
                body.append(text(x + bar_width / 2, bottom + 25 + index * 17, part, size=14))
        body.append(text(left - 56, top + panel_height / 2, "Saturated windows (%)", size=15, rotate=-90))

    out_base = run_dir / "comparison_to_manuscript" / "published_vs_eigenvalue_weighted_saturated_windows"
    convert_svg(svg_document(width, height, body), out_base, width)


def main() -> None:
    args = parse_args()
    run_dir = args.run_dir.resolve()
    rows = read_summary(run_dir)
    plot_curves(rows, run_dir)
    plot_saturated_windows(rows, run_dir)


if __name__ == "__main__":
    main()
