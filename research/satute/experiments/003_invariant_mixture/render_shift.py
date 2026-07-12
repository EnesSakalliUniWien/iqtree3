#!/usr/bin/env python3
"""Render direct paired shifts from weighted_invariant_simulations.py output."""

import argparse
import csv
import math
import shutil
import subprocess
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]


MODEL_STYLE = {
    "+I+G4": ("#0A3FBA", ""),
    "+I+R4": ("#B6497D", "7,4"),
}


def read_rows(path, text_fields):
    with path.open("r", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        for key, value in list(row.items()):
            if key not in text_fields:
                row[key] = float(value)
    return rows


def scale_log(value, lower, upper, start, length):
    return start + length * (math.log10(value) - math.log10(lower)) / (
        math.log10(upper) - math.log10(lower)
    )


def scale_y(value, lower, upper, start, length):
    return start + length * (1.0 - (value - lower) / (upper - lower))


def padded_range(rows, keys, include_zero=True):
    values = [row[key] for row in rows for key in keys]
    if include_zero:
        values.append(0.0)
    lower, upper = min(values), max(values)
    span = max(upper - lower, 1e-6)
    return lower - 0.10 * span, upper + 0.10 * span


def axes(svg, left, top, width, height, title, subtitle, ymin, ymax, ylabel):
    svg.append(
        f'<text x="{left}" y="{top-28}" font-size="15" font-weight="700" '
        f'fill="#202124">{title}</text>'
    )
    svg.append(
        f'<text x="{left}" y="{top-10}" font-size="10.5" fill="#5f6368">{subtitle}</text>'
    )
    svg.append(
        f'<rect x="{left}" y="{top}" width="{width}" height="{height}" '
        'fill="#ffffff" stroke="#3c4043"/>'
    )
    for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
        value = ymin + fraction * (ymax - ymin)
        y = scale_y(value, ymin, ymax, top, height)
        svg.append(
            f'<line x1="{left}" x2="{left+width}" y1="{y:.2f}" y2="{y:.2f}" '
            'stroke="#eceff1"/>'
        )
        svg.append(
            f'<text x="{left-7}" y="{y+4:.2f}" text-anchor="end" font-size="9" '
            f'fill="#5f6368">{value:.3g}</text>'
        )
    if ymin <= 0.0 <= ymax:
        y = scale_y(0.0, ymin, ymax, top, height)
        svg.append(
            f'<line x1="{left}" x2="{left+width}" y1="{y:.2f}" y2="{y:.2f}" '
            'stroke="#3c4043" stroke-dasharray="4,4"/>'
        )
    for tick in (100, 500, 1000, 5000, 10000):
        x = scale_log(tick, 100, 10000, left, width)
        svg.append(
            f'<text x="{x:.2f}" y="{top+height+16}" text-anchor="middle" '
            f'font-size="8.5" fill="#5f6368">{tick:g}</text>'
        )
    svg.append(
        f'<text x="{left+width/2}" y="{top+height+36}" text-anchor="middle" '
        'font-size="10" fill="#3c4043">alignment length, n (log scale)</text>'
    )
    svg.append(
        f'<text x="{left-55}" y="{top+height/2}" text-anchor="middle" '
        f'font-size="10" fill="#3c4043" transform="rotate(-90 {left-55} '
        f'{top+height/2})">{ylabel}</text>'
    )


def line_with_band(svg, rows, left, top, width, height, ymin, ymax, center, low, high, color, dash):
    rows = sorted(rows, key=lambda row: row["nsites"])
    upper = [
        (
            scale_log(row["nsites"], 100, 10000, left, width),
            scale_y(row[high], ymin, ymax, top, height),
        )
        for row in rows
    ]
    lower = [
        (
            scale_log(row["nsites"], 100, 10000, left, width),
            scale_y(row[low], ymin, ymax, top, height),
        )
        for row in reversed(rows)
    ]
    points = " ".join(f"{x:.2f},{y:.2f}" for x, y in upper + lower)
    svg.append(f'<polygon points="{points}" fill="{color}" fill-opacity="0.11"/>')
    centers = [
        (
            scale_log(row["nsites"], 100, 10000, left, width),
            scale_y(row[center], ymin, ymax, top, height),
        )
        for row in rows
    ]
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    svg.append(
        '<polyline points="'
        + " ".join(f"{x:.2f},{y:.2f}" for x, y in centers)
        + f'" fill="none" stroke="{color}" stroke-width="2.3"{dash_attr}/>'
    )
    for x, y in centers:
        svg.append(
            f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2.8" fill="#fff" '
            f'stroke="{color}" stroke-width="1.5"/>'
        )


def render(paired, output_prefix, null_reps, alternative_reps, png_scale=3.0):
    width, height = 1240, 860
    panel_w, panel_h = 480, 250
    panels = ((90, 160), (675, 160), (90, 550), (675, 550))
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="620" y="38" text-anchor="middle" font-size="24" font-weight="700" fill="#202124">How much the refined +I statistic moves</text>',
        f'<text x="620" y="64" text-anchor="middle" font-size="12" fill="#5f6368">Paired exact-pattern rerun; {null_reps:,} null and {alternative_reps:,} alternative replicates per cell</text>',
    ]

    selections = [
        (
            "A  Null shift in standardized statistic",
            "t = 10,000; refined minus legacy on the same null alignments",
            [row for row in paired if row["sample"] == "null_heldout" and row["rule"] == "nominal" and row["branch_length"] == 10000],
            "mean_delta_z",
            "mean_delta_z_ci_low",
            "mean_delta_z_ci_high",
            "paired mean ΔZ",
        ),
        (
            "B  Nominal false-positive shift",
            "t = 10,000; refined minus legacy rejection fraction",
            [row for row in paired if row["sample"] == "null_heldout" and row["rule"] == "nominal" and row["branch_length"] == 10000],
            "refined_minus_legacy",
            "paired_bootstrap_ci_low",
            "paired_bootstrap_ci_high",
            "Δ false-positive rate",
        ),
        (
            "C  Calibrated decision-margin shift",
            "+I+G4 at t = 50; +I+R4 at t = 20",
            [row for row in paired if row["sample"] == "alternative" and row["rule"] == "calibrated" and ((row["model"] == "+I+G4" and row["branch_length"] == 50) or (row["model"] == "+I+R4" and row["branch_length"] == 20))],
            "mean_delta_decision_margin",
            "mean_delta_decision_margin_ci_low",
            "mean_delta_decision_margin_ci_high",
            "Δ(Z − calibrated threshold)",
        ),
        (
            "D  Calibrated-power shift",
            "+I+G4 at t = 50; +I+R4 at t = 20",
            [row for row in paired if row["sample"] == "alternative" and row["rule"] == "calibrated" and ((row["model"] == "+I+G4" and row["branch_length"] == 50) or (row["model"] == "+I+R4" and row["branch_length"] == 20))],
            "refined_minus_legacy",
            "paired_bootstrap_ci_low",
            "paired_bootstrap_ci_high",
            "Δ calibrated power",
        ),
    ]

    for panel, selection in zip(panels, selections):
        left, top = panel
        title, subtitle, rows, center, low, high, ylabel = selection
        ymin, ymax = padded_range(rows, (low, high))
        axes(svg, left, top, panel_w, panel_h, title, subtitle, ymin, ymax, ylabel)
        for model, (color, dash) in MODEL_STYLE.items():
            model_rows = [row for row in rows if row["model"] == model]
            line_with_band(
                svg,
                model_rows,
                left,
                top,
                panel_w,
                panel_h,
                ymin,
                ymax,
                center,
                low,
                high,
                color,
                dash,
            )

    for index, (model, (color, dash)) in enumerate(MODEL_STYLE.items()):
        x = 440 + index * 260
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        svg.append(
            f'<line x1="{x}" x2="{x+40}" y1="92" y2="92" stroke="{color}" '
            f'stroke-width="2.5"{dash_attr}/>'
        )
        svg.append(
            f'<text x="{x+50}" y="96" font-size="11" fill="#3c4043">{model}</text>'
        )
    svg.append(
        '<text x="90" y="842" font-size="10" fill="#5f6368">Negative null shifts remove legacy bias; positive decision-margin and calibrated-power shifts favor the refined statistic.</text>'
    )
    svg.append("</svg>")

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    svg_path = output_prefix.with_suffix(".svg")
    svg_path.write_text("\n".join(svg) + "\n", encoding="utf-8")
    rsvg = shutil.which("rsvg-convert")
    if rsvg:
        subprocess.run(
            [rsvg, "-f", "png", "-z", str(png_scale), "-o", str(output_prefix.with_suffix(".png")), str(svg_path)],
            check=True,
        )
        subprocess.run(
            [rsvg, "-f", "pdf", "-o", str(output_prefix.with_suffix(".pdf")), str(svg_path)],
            check=True,
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-dir",
        default=str(PROJECT_ROOT / "artifacts/local/weighted-invariant-shift"),
    )
    parser.add_argument(
        "--output-prefix",
        default=str(PROJECT_ROOT / "artifacts/local/weighted-invariant-shift/figures/direct-shift"),
    )
    parser.add_argument("--png-scale", type=float, default=3.0)
    args = parser.parse_args()
    input_dir = Path(args.input_dir)
    paired = read_rows(
        input_dir / "weighted_invariant_paired.tsv",
        {"model", "rule", "sample"},
    )
    null_reps = int(
        max(row["replicates"] for row in paired if row["sample"] == "null_heldout") * 2
    )
    alternative_reps = int(
        max(row["replicates"] for row in paired if row["sample"] == "alternative")
    )
    render(paired, Path(args.output_prefix), null_reps, alternative_reps, args.png_scale)
    print(Path(args.output_prefix).with_suffix(".svg"))


if __name__ == "__main__":
    main()
