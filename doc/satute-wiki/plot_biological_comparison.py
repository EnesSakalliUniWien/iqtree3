#!/usr/bin/env python3
"""Regenerate biological sliding-window comparison figures from curated tables."""

from __future__ import annotations

import argparse
import csv
import html
import math
import shutil
import subprocess
import tempfile
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DEFAULT_RUN_DIR = (
    SCRIPT_DIR
    / "results"
    / "biological"
    / "enhanced_sliding_window_rerun_20260619_150104"
)

BRANCH_ORDER = ("eukaryota", "yeast")
FORMULA_ORDER = ("original_satute", "dominant", "eigenvalue_weighted")
LABELS = {
    "original_satute": "published SatuTe",
    "dominant": "dominant",
    "eigenvalue_weighted": "eigenvalue-weighted",
}
COLORS = {
    "original_satute": "#222222",
    "dominant": "#4c78a8",
    "eigenvalue_weighted": "#54a24b",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot curated biological SatuTe sliding-window comparisons."
    )
    parser.add_argument(
        "--run-dir",
        type=Path,
        default=DEFAULT_RUN_DIR,
        help="Curated biological rerun directory.",
    )
    return parser.parse_args()


def resolve_data_path(path_text: str, run_dir: Path) -> Path:
    path = Path(path_text)
    candidates: list[Path] = []
    if path.is_absolute():
        candidates.append(path)
    else:
        candidates.extend((Path.cwd() / path, REPO_ROOT / path, run_dir / path))
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(path_text)


def read_summary(run_dir: Path) -> list[dict[str, str]]:
    summary_path = run_dir / "comparison_to_manuscript" / "original_vs_enhanced_summary.tsv"
    with summary_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    expected = {(formula, branch) for formula in FORMULA_ORDER for branch in BRANCH_ORDER}
    observed = {(row["formula"], row["branch_short"]) for row in rows}
    missing = sorted(expected - observed)
    extra = sorted(observed - expected)
    if missing or extra:
        raise ValueError(f"Unexpected summary formulas or branches: missing={missing}, extra={extra}")
    return rows


def read_window_scores(csv_path: Path, expected_windows: int) -> list[float]:
    scores: list[float] = []
    with csv_path.open() as handle:
        for raw_line in handle:
            value = raw_line.strip().strip('"')
            if not value:
                continue
            try:
                scores.append(float(value))
            except ValueError:
                continue

    if len(scores) != expected_windows:
        raise ValueError(
            f"{csv_path} contains {len(scores)} scores, expected {expected_windows}"
        )
    return scores


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


def text(
    x: float,
    y: float,
    value: str,
    *,
    size: int = 18,
    anchor: str = "middle",
    weight: str = "400",
    rotate: float | None = None,
) -> str:
    transform = f' transform="rotate({rotate:.1f} {x:.2f} {y:.2f})"' if rotate is not None else ""
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" text-anchor="{anchor}" '
        f'font-family="Arial, Helvetica, sans-serif" font-size="{size}" '
        f'font-weight="{weight}" fill="#111111"{transform}>{esc(value)}</text>'
    )


def line(x1: float, y1: float, x2: float, y2: float, *, color: str, width: float = 1.0, dash: str = "") -> str:
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    return (
        f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
        f'stroke="{color}" stroke-width="{width:.2f}"{dash_attr}/>'
    )


def rect(x: float, y: float, width: float, height: float, *, fill: str, stroke: str = "none") -> str:
    return (
        f'<rect x="{x:.2f}" y="{y:.2f}" width="{width:.2f}" height="{height:.2f}" '
        f'fill="{fill}" stroke="{stroke}"/>'
    )


def polyline(points: list[tuple[float, float]], *, color: str, width: float) -> str:
    coords = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
    return (
        f'<polyline points="{coords}" fill="none" stroke="{color}" '
        f'stroke-width="{width:.2f}" stroke-linejoin="round" stroke-linecap="round"/>'
    )


def svg_document(width: int, height: int, body: list[str]) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">\n'
        '<rect width="100%" height="100%" fill="#ffffff"/>\n'
        + "\n".join(body)
        + "\n</svg>\n"
    )


def nice_ticks(y_min: float, y_max: float, step: float) -> list[float]:
    start = step * math.floor(y_min / step)
    end = step * math.ceil(y_max / step)
    ticks = []
    value = start
    while value <= end + step / 2:
        ticks.append(round(value, 6))
        value += step
    return ticks


def convert_svg(svg: str, out_base: Path, width: int) -> None:
    converter = shutil.which("rsvg-convert")
    if converter is None:
        raise RuntimeError("rsvg-convert is required to write PDF and PNG figures")

    with tempfile.NamedTemporaryFile("w", suffix=".svg", delete=False) as handle:
        handle.write(svg)
        svg_path = Path(handle.name)
    try:
        subprocess.run(
            [converter, "-f", "pdf", "-o", str(out_base.with_suffix(".pdf")), str(svg_path)],
            check=True,
        )
        subprocess.run(
            [
                converter,
                "-f",
                "png",
                "-w",
                str(width),
                "-o",
                str(out_base.with_suffix(".png")),
                str(svg_path),
            ],
            check=True,
        )
    finally:
        svg_path.unlink(missing_ok=True)


def plot_curves(rows: list[dict[str, str]], run_dir: Path) -> None:
    by_key = {(row["formula"], row["branch_short"]): row for row in rows}
    series: dict[tuple[str, str], list[float]] = {}
    for branch in BRANCH_ORDER:
        for formula in FORMULA_ORDER:
            row = by_key[(formula, branch)]
            series[(formula, branch)] = read_window_scores(
                resolve_data_path(row["csv"], run_dir),
                int(row["windows"]),
            )

    values = [score for scores in series.values() for score in scores]
    y_min = min(-1.4, min(values) - 0.4)
    y_max = max(values) + 0.8
    ticks = nice_ticks(y_min, y_max, 2.0)

    width, height = 1650, 900
    plot_left, plot_right = 95, width - 35
    panel_top = (115, 515)
    panel_height = 295
    plot_width = plot_right - plot_left
    body: list[str] = []

    legend_y = 38
    legend_x = 560
    for index, formula in enumerate(FORMULA_ORDER):
        x = legend_x + index * 250
        body.append(line(x, legend_y - 6, x + 58, legend_y - 6, color=COLORS[formula], width=3.0))
        body.append(text(x + 72, legend_y, LABELS[formula], size=20, anchor="start"))

    for branch_index, branch in enumerate(BRANCH_ORDER):
        top = panel_top[branch_index]
        bottom = top + panel_height
        body.append(text(width / 2, top - 24, branch, size=28, weight="700"))
        body.append(rect(plot_left, top, plot_width, panel_height, fill="#ffffff", stroke="#111111"))

        def sx(index: int, n_values: int) -> float:
            return plot_left + plot_width * index / max(1, n_values - 1)

        def sy(value: float) -> float:
            return bottom - panel_height * (value - y_min) / (y_max - y_min)

        for tick in ticks:
            if y_min <= tick <= y_max:
                y = sy(tick)
                body.append(line(plot_left, y, plot_right, y, color="#d7d7d7", width=0.8))
                body.append(text(plot_left - 14, y + 6, f"{tick:g}", size=16, anchor="end"))

        n_windows = int(by_key[(FORMULA_ORDER[0], branch)]["windows"])
        for tick in range(0, n_windows, 250):
            x = sx(tick, n_windows)
            body.append(line(x, top, x, bottom, color="#eeeeee", width=0.8))
            if branch_index == len(BRANCH_ORDER) - 1:
                body.append(text(x, bottom + 24, str(tick), size=15))

        threshold = float(by_key[(FORMULA_ORDER[0], branch)]["z_alpha"])
        body.append(line(plot_left, sy(threshold), plot_right, sy(threshold), color="#999999", width=1.5, dash="8 6"))

        for formula in FORMULA_ORDER:
            scores = series[(formula, branch)]
            points = [(sx(i, len(scores)), sy(value)) for i, value in enumerate(scores)]
            body.append(polyline(points, color=COLORS[formula], width=1.35))

        body.append(text(32, top + panel_height / 2, "Window z-score", size=18, rotate=-90))

    body.append(text(width / 2, height - 16, "Window index", size=20))
    svg = svg_document(width, height, body)
    out_base = run_dir / "comparison_to_manuscript" / "original_vs_enhanced_sliding_window_curves"
    convert_svg(svg, out_base, width)


def plot_saturated_windows(rows: list[dict[str, str]], run_dir: Path) -> None:
    by_key = {(row["formula"], row["branch_short"]): row for row in rows}
    width, height = 1100, 510
    plot_top, plot_bottom = 80, 420
    panel_width = 430
    panel_lefts = (95, 620)
    y_max = 52.0
    body: list[str] = []

    for branch, left in zip(BRANCH_ORDER, panel_lefts, strict=True):
        body.append(text(left + panel_width / 2, 42, branch, size=28, weight="700"))
        body.append(rect(left, plot_top, panel_width, plot_bottom - plot_top, fill="#ffffff", stroke="#111111"))

        for tick in range(0, 51, 10):
            y = plot_bottom - (plot_bottom - plot_top) * tick / y_max
            body.append(line(left, y, left + panel_width, y, color="#d7d7d7", width=0.9))
            body.append(text(left - 12, y + 6, str(tick), size=16, anchor="end"))

        bar_width = 82
        gap = 40
        total = len(FORMULA_ORDER) * bar_width + (len(FORMULA_ORDER) - 1) * gap
        start = left + (panel_width - total) / 2
        for index, formula in enumerate(FORMULA_ORDER):
            value = 100.0 * float(by_key[(formula, branch)]["saturated_fraction"])
            bar_height = (plot_bottom - plot_top) * value / y_max
            x = start + index * (bar_width + gap)
            y = plot_bottom - bar_height
            body.append(rect(x, y, bar_width, bar_height, fill=COLORS[formula]))
            label_parts = LABELS[formula].split(" ")
            for part_index, part in enumerate(label_parts):
                body.append(text(x + bar_width / 2, plot_bottom + 28 + part_index * 18, part, size=15))

        body.append(text(left - 64, (plot_top + plot_bottom) / 2, "Saturated windows (%)", size=18, rotate=-90))

    svg = svg_document(width, height, body)
    out_base = run_dir / "comparison_to_manuscript" / "original_vs_enhanced_saturated_windows"
    convert_svg(svg, out_base, width)


def main() -> None:
    args = parse_args()
    run_dir = args.run_dir.resolve()
    rows = read_summary(run_dir)
    plot_curves(rows, run_dir)
    plot_saturated_windows(rows, run_dir)


if __name__ == "__main__":
    main()
