#!/usr/bin/env python3
"""Render the original-style six-panel model-misspecification comparison."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path
from xml.sax.saxutils import escape


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from satute_analysis.publication_figures import (  # noqa: E402
    MISSPECIFICATION_CONTRACT,
    stable_case_summary,
    stable_figure_root,
)
from tools.figures.export import export_svg_pdf, export_svg_png  # noqa: E402


OUTDIR = stable_figure_root(MISSPECIFICATION_CONTRACT.id)
TREE_LABELS = {
    "five_external": "five-taxon external branch",
    "sixteen_internal": "16-taxon internal branch",
}
MODEL_LABELS = {"JC": "JC evaluation", "K2P": "K2P evaluation", "F81": "F81 evaluation"}
SERIES_COLORS = {
    "correct_reference": "#E69F00",
    "misspecified_true_tree": "#D55E72",
    "misspecified_ml_lengths": "#8E2C68",
    "misspecified_ml_tree_bonferroni": "#007A3D",
}
SITE_DASHES = {100: "10 6", 250: "4 5", 1000: ""}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--summary",
        type=Path,
        default=stable_case_summary("misspec_JC"),
        help="Stable schema-v2 misspecification summary.",
    )
    parser.add_argument("--output-dir", type=Path, default=OUTDIR)
    parser.add_argument("--expected-reps", type=int, default=1000)
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "tree_case",
        "simulation_model",
        "evaluation_model",
        "nsites",
        "branch_length",
        "scenario",
        "decision_rule",
        "formula",
        "evaluated",
        "informative",
        "missing_split",
        "fraction_informative",
    }
    missing = required - set(rows[0] if rows else ())
    if missing:
        raise SystemExit(f"Summary is missing columns: {sorted(missing)}")
    return rows


def selected_rows(
    rows: list[dict[str, str]], expected_reps: int
) -> tuple[list[dict[str, str]], str]:
    simulation_models = {row["simulation_model"] for row in rows}
    if len(simulation_models) != 1:
        raise SystemExit("Misspecification summary must contain one simulation model")
    simulation_model = next(iter(simulation_models))
    output: list[dict[str, str]] = []
    for evaluation_model in MISSPECIFICATION_CONTRACT.evaluation_models:
        for series in MISSPECIFICATION_CONTRACT.series:
            target_evaluation = (
                simulation_model if series.reference_model else evaluation_model
            )
            selected = [
                row
                for row in rows
                if row["simulation_model"] == simulation_model
                and row["evaluation_model"] == target_evaluation
                and row["formula"] == "dominant"
                and row["scenario"] == series.scenario
                and row["decision_rule"] == series.decision_rule
                and row["tree_case"] in MISSPECIFICATION_CONTRACT.tree_cases
            ]
            if not selected:
                raise SystemExit(
                    f"No rows for {evaluation_model}/{series.scenario}/"
                    f"{series.decision_rule}"
                )
            for row in selected:
                total = int(row["evaluated"]) + int(row["missing_split"])
                if total != expected_reps:
                    raise SystemExit(
                        f"Incomplete cell for {evaluation_model}/{series.id}: {total}"
                    )
                output.append(
                    {
                        **row,
                        "panel_evaluation_model": evaluation_model,
                        "series_id": series.id,
                        "series_label": series.label,
                    }
                )
    return output, simulation_model


def scale(value: float, low: float, high: float, start: float, length: float) -> float:
    if math.isclose(low, high):
        return start + length / 2
    return start + (value - low) / (high - low) * length


def path_data(points: list[tuple[float, float]]) -> str:
    return " ".join(
        [f"M {points[0][0]:.2f} {points[0][1]:.2f}"]
        + [f"L {x:.2f} {y:.2f}" for x, y in points[1:]]
    )


def render(rows: list[dict[str, str]], output_dir: Path) -> None:
    width, height = 1500, 1710
    panel_width, panel_height = 570, 280
    x_origins = (145, 815)
    y_origins = (150, 565, 980)
    branch_values = sorted({float(row["branch_length"]) for row in rows})
    site_values = sorted({int(row["nsites"]) for row in rows})
    xmin, xmax = min(branch_values), max(branch_values)
    body = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#202124}</style>',
        f'<text x="{width / 2}" y="42" text-anchor="middle" font-size="28" font-weight="700">SatuTe under substitution-model misspecification</text>',
        f'<text x="{width / 2}" y="74" text-anchor="middle" font-size="16" fill="#4a4a4a">GTR+F simulations; 1,000 replicates per point; published SatuTe statistic</text>',
    ]
    grouped: dict[tuple[str, str, str, int], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[
            (
                row["panel_evaluation_model"],
                row["tree_case"],
                row["series_id"],
                int(row["nsites"]),
            )
        ].append(row)

    letter = 0
    for model_index, evaluation_model in enumerate(
        MISSPECIFICATION_CONTRACT.evaluation_models
    ):
        for tree_index, tree_case in enumerate(MISSPECIFICATION_CONTRACT.tree_cases):
            left, top = x_origins[tree_index], y_origins[model_index]
            body.append(
                f'<rect x="{left}" y="{top}" width="{panel_width}" height="{panel_height}" '
                'fill="#ffffff" stroke="#3f3f3f" stroke-width="1.2"/>'
            )
            band_top = top + panel_height * 0.95
            body.append(
                f'<rect x="{left}" y="{band_top:.2f}" width="{panel_width}" '
                f'height="{panel_height * 0.05:.2f}" fill="#e0e0e0"/>'
            )
            body.append(
                f'<text x="{left - 42}" y="{top + 21}" font-size="25" font-weight="700">'
                f'{chr(ord("a") + letter)}</text>'
            )
            letter += 1
            body.append(
                f'<text x="{left + panel_width / 2}" y="{top - 16}" text-anchor="middle" '
                f'font-size="18" font-weight="700">{escape(MODEL_LABELS[evaluation_model])} — '
                f'{escape(TREE_LABELS[tree_case])}</text>'
            )
            for value in (0, 25, 50, 75, 100):
                y = top + panel_height * (1 - value / 100)
                body.append(
                    f'<line x1="{left}" x2="{left + panel_width}" y1="{y:.2f}" y2="{y:.2f}" '
                    'stroke="#e7e7e7" stroke-width="1"/>'
                )
                if tree_index == 0:
                    body.append(
                        f'<text x="{left - 12}" y="{y + 5:.2f}" text-anchor="end" '
                        f'font-size="14">{value}</text>'
                    )
            for branch in branch_values:
                x = scale(branch, xmin, xmax, left, panel_width)
                body.append(
                    f'<line x1="{x:.2f}" x2="{x:.2f}" y1="{top}" y2="{top + panel_height}" '
                    'stroke="#f0f0f0" stroke-width="0.8"/>'
                )
                if branch in (4, 6, 8, 10, 12):
                    body.append(
                        f'<text x="{x:.2f}" y="{top + panel_height + 24}" '
                        f'text-anchor="middle" font-size="14">{branch:g}</text>'
                    )
            for series in MISSPECIFICATION_CONTRACT.series:
                for nsites in site_values:
                    selected = sorted(
                        grouped[(evaluation_model, tree_case, series.id, nsites)],
                        key=lambda row: float(row["branch_length"]),
                    )
                    points = [
                        (
                            scale(float(row["branch_length"]), xmin, xmax, left, panel_width),
                            top
                            + panel_height
                            * (1 - 100 * float(row["fraction_informative"]) / 100),
                        )
                        for row in selected
                    ]
                    if not points:
                        continue
                    dash = SITE_DASHES.get(nsites, "")
                    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
                    color = SERIES_COLORS[series.id]
                    body.append(
                        f'<path d="{path_data(points)}" fill="none" stroke="{color}" '
                        f'stroke-width="2.1"{dash_attr}/>'
                    )
                    for x, y in points:
                        body.append(
                            f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2.6" '
                            f'fill="{color}" stroke="#ffffff" stroke-width="0.7"/>'
                        )
            if tree_index == 0:
                body.append(
                    f'<text x="{left - 72}" y="{top + panel_height / 2}" text-anchor="middle" '
                    f'font-size="16" transform="rotate(-90 {left - 72} {top + panel_height / 2})">'
                    "phylogenetically informative alignments (%)</text>"
                )
            body.append(
                f'<text x="{left + panel_width / 2}" y="{top + panel_height + 50}" '
                'text-anchor="middle" font-size="16">target branch length</text>'
            )

    legend_y = 1385
    body.append(
        f'<text x="145" y="{legend_y}" font-size="18" font-weight="700">Analysis series</text>'
    )
    for index, series in enumerate(MISSPECIFICATION_CONTRACT.series):
        x = 145 + (index % 2) * 670
        y = legend_y + 34 + (index // 2) * 40
        color = SERIES_COLORS[series.id]
        body.append(
            f'<line x1="{x}" x2="{x + 48}" y1="{y}" y2="{y}" stroke="{color}" stroke-width="3"/>'
        )
        body.append(
            f'<text x="{x + 60}" y="{y + 6}" font-size="16">{escape(series.label)}</text>'
        )
    site_y = 1515
    body.append(
        f'<text x="145" y="{site_y}" font-size="18" font-weight="700">Alignment length</text>'
    )
    for index, nsites in enumerate(site_values):
        x = 365 + index * 285
        dash = SITE_DASHES.get(nsites, "")
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        body.append(
            f'<line x1="{x}" x2="{x + 55}" y1="{site_y - 5}" y2="{site_y - 5}" '
            f'stroke="#222222" stroke-width="2.6"{dash_attr}/>'
        )
        body.append(
            f'<text x="{x + 67}" y="{site_y + 1}" font-size="16">{nsites:,} sites</text>'
        )
    body.append(
        '<text x="145" y="1640" font-size="14" fill="#555555">'
        "Gray band marks the nominal 0–5% rejection region. "
        "Bonferroni is used only for the inferred-ML-tree series.</text>"
    )
    body.append("</svg>")

    output_dir.mkdir(parents=True, exist_ok=True)
    svg_path = output_dir / "figure_model_misspecification.svg"
    pdf_path = output_dir / "figure_model_misspecification.pdf"
    png_path = output_dir / "figure_model_misspecification.png"
    svg_path.write_text("\n".join(body) + "\n", encoding="utf-8")
    export_svg_pdf(svg_path, pdf_path)
    export_svg_png(svg_path, png_path, width=3000)


def write_selected(rows: list[dict[str, str]], output_dir: Path) -> None:
    path = output_dir / "figure_model_misspecification_data.tsv"
    fields = [
        "panel_evaluation_model",
        "tree_case",
        "series_id",
        "series_label",
        "simulation_model",
        "evaluation_model",
        "scenario",
        "decision_rule",
        "formula",
        "nsites",
        "branch_length",
        "evaluated",
        "informative",
        "missing_split",
        "fraction_informative",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows({field: row[field] for field in fields} for row in rows)


def main() -> None:
    args = parse_args()
    rows, _simulation_model = selected_rows(read_rows(args.summary), args.expected_reps)
    render(rows, args.output_dir / "figures")
    (args.output_dir / "data").mkdir(parents=True, exist_ok=True)
    write_selected(rows, args.output_dir / "data")
    print(f"Rendered {len(rows)} reviewed Figure 5 summary cells in {args.output_dir}")


if __name__ == "__main__":
    main()
