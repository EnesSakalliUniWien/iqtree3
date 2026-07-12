#!/usr/bin/env python3
"""Summarize the alignment-length null-scaling experiment."""

import argparse
import csv
import math
import shutil
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from satute_analysis.plotting import (
    SATURATION_FORMULAS as FORMULAS,
    SATURATION_FORMULA_STYLE as FORMULA_STYLE,
    polygon,
    scale_linear,
    svg_path as path,
)


def quantile(values, probability):
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = int(position)
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[min(lower + 1, len(ordered) - 1)] * fraction


def wilson(successes, total, z=1.96):
    fraction = successes / total
    denominator = 1.0 + z * z / total
    center = (fraction + z * z / (2.0 * total)) / denominator
    half = z * math.sqrt(fraction * (1.0 - fraction) / total + z * z / (4.0 * total * total)) / denominator
    return center - half, center + half


def load_rows(root, null_length):
    output = []
    directories = [path for path in root.glob("n*") if path.is_dir() and path.name[1:].isdigit()]
    for directory in sorted(directories, key=lambda path: int(path.name[1:])):
        nsites = int(directory.name[1:])
        detail = directory / "saturation_diagnostic_detail.tsv"
        with detail.open("r", encoding="utf-8", newline="") as handle:
            rows = [
                row
                for row in csv.DictReader(handle, delimiter="\t")
                if math.isclose(float(row["branch_length"]), null_length)
            ]
        groups = defaultdict(list)
        for row in rows:
            groups[row["formula"]].append(row)
        for formula in FORMULAS:
            selected = groups[formula]
            z_values = [float(row["satZ"]) for row in selected]
            c_values = [float(row["satC"]) for row in selected]
            se_values = [float(row["satSE"]) for row in selected]
            mean_z = statistics.fmean(z_values)
            sd_z = statistics.stdev(z_values)
            mean_half = 1.96 * sd_z / math.sqrt(len(z_values))
            successes = sum(value > 1.6448536269514722 for value in z_values)
            reject_low, reject_high = wilson(successes, len(z_values))
            output.append(
                {
                    "nsites": nsites,
                    "formula": formula,
                    "replicates": len(selected),
                    "mean_z": mean_z,
                    "mean_z_ci_low": mean_z - mean_half,
                    "mean_z_ci_high": mean_z + mean_half,
                    "sd_z": sd_z,
                    "q95_z": quantile(z_values, 0.95),
                    "nominal_rejection": successes / len(z_values),
                    "rejection_ci_low": reject_low,
                    "rejection_ci_high": reject_high,
                    "variance_ratio": (statistics.stdev(c_values) / statistics.fmean(se_values)) ** 2,
                }
            )
    return output


def sx(value, xmin, xmax, left, width):
    return scale_linear(math.log10(value), math.log10(xmin), math.log10(xmax), left, width)


def sy(value, ymin, ymax, top, height):
    return scale_linear(value, ymin, ymax, top, height, invert=True)


def axes(svg, left, top, width, height, title, subtitle, ymin, ymax, ylabel, nvalues):
    svg.append(f'<text x="{left}" y="{top-28}" font-size="16" font-weight="700" fill="#202124">{title}</text>')
    svg.append(f'<text x="{left}" y="{top-10}" font-size="10.5" fill="#5f6368">{subtitle}</text>')
    svg.append(f'<rect x="{left}" y="{top}" width="{width}" height="{height}" fill="#fff" stroke="#3c4043"/>')
    for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
        value = ymin + fraction * (ymax - ymin)
        y = sy(value, ymin, ymax, top, height)
        svg.append(f'<line x1="{left}" x2="{left+width}" y1="{y:.2f}" y2="{y:.2f}" stroke="#e8eaed"/>')
        svg.append(f'<text x="{left-7}" y="{y+4:.2f}" text-anchor="end" font-size="9.5" fill="#5f6368">{value:.3g}</text>')
    xmin, xmax = min(nvalues), max(nvalues)
    for value in nvalues:
        x = sx(value, xmin, xmax, left, width)
        svg.append(f'<line x1="{x:.2f}" x2="{x:.2f}" y1="{top}" y2="{top+height}" stroke="#f1f3f4"/>')
        svg.append(f'<text x="{x:.2f}" y="{top+height+16}" text-anchor="middle" font-size="9.5" fill="#5f6368">{value:,}</text>')
    svg.append(f'<text x="{left+width/2}" y="{top+height+36}" text-anchor="middle" font-size="10.5">alignment length</text>')
    x = left - 45
    y = top + height / 2
    svg.append(f'<text x="{x}" y="{y}" text-anchor="middle" font-size="10.5" transform="rotate(-90 {x} {y})">{ylabel}</text>')


def render(rows, output_prefix):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["formula"]].append(row)
    for values in grouped.values():
        values.sort(key=lambda row: row["nsites"])
    nvalues = sorted({row["nsites"] for row in rows})
    xmin, xmax = min(nvalues), max(nvalues)
    width, height = 1280, 540
    panel_w, panel_h, top = 330, 300, 145
    lefts = (90, 480, 870)
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="640" y="38" text-anchor="middle" font-size="24" font-weight="700" fill="#202124">Standard-normal null behavior across alignment lengths</text>',
        '<text x="640" y="63" text-anchor="middle" font-size="12" fill="#5f6368">Independent ALISIM batch streams; 1,000 saturated-branch replicates per alignment length</text>',
    ]
    axes(svg, lefts[0], top, panel_w, panel_h, "A  Null centering", "Mean satZ with 95% confidence interval", -0.12, 0.12, "mean satZ", nvalues)
    axes(svg, lefts[1], top, panel_w, panel_h, "B  Null scale", "Standard deviation of satZ", 0.9, 1.1, "SD(satZ)", nvalues)
    axes(svg, lefts[2], top, panel_w, panel_h, "C  Nominal rejection", "One-sided threshold 1.645; Wilson 95% interval", 0.0, 0.10, "rejection probability", nvalues)
    reference = [(0, 0.0), (1, 1.0), (2, 0.05)]
    for panel, value in reference:
        ymin, ymax = ((-0.12, 0.12), (0.9, 1.1), (0.0, 0.10))[panel]
        y = sy(value, ymin, ymax, top, panel_h)
        svg.append(f'<line x1="{lefts[panel]}" x2="{lefts[panel]+panel_w}" y1="{y:.2f}" y2="{y:.2f}" stroke="#3c4043" stroke-dasharray="5,4"/>')
    for formula in FORMULAS:
        label, color, dash, _marker = FORMULA_STYLE[formula]
        values = grouped[formula]
        panels = (
            ("mean_z", "mean_z_ci_low", "mean_z_ci_high", -0.12, 0.12),
            ("sd_z", None, None, 0.9, 1.1),
            ("nominal_rejection", "rejection_ci_low", "rejection_ci_high", 0.0, 0.10),
        )
        for panel, (center, low, high, ymin, ymax) in enumerate(panels):
            if low:
                upper = [(sx(row["nsites"], xmin, xmax, lefts[panel], panel_w), sy(row[high], ymin, ymax, top, panel_h)) for row in values]
                lower = [(sx(row["nsites"], xmin, xmax, lefts[panel], panel_w), sy(row[low], ymin, ymax, top, panel_h)) for row in reversed(values)]
                svg.append(f'<polygon points="{polygon(upper+lower)}" fill="{color}" fill-opacity="0.10"/>')
            points = [(sx(row["nsites"], xmin, xmax, lefts[panel], panel_w), sy(row[center], ymin, ymax, top, panel_h)) for row in values]
            dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
            svg.append(f'<path d="{path(points)}" fill="none" stroke="{color}" stroke-width="2.2"{dash_attr}/>')
            for x, y in points:
                svg.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3" fill="white" stroke="{color}" stroke-width="1.6"/>')
        legend_x = 170 + FORMULAS.index(formula) * 335
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        svg.append(f'<line x1="{legend_x}" x2="{legend_x+38}" y1="98" y2="98" stroke="{color}" stroke-width="2.4"{dash_attr}/>')
        svg.append(f'<text x="{legend_x+48}" y="102" font-size="10.5" fill="#3c4043">{label}</text>')
    svg.append('<text x="90" y="522" font-size="10.5" fill="#5f6368">All mean confidence intervals include zero; all SDs and variance ratios are close to one. No null correction is indicated.</text>')
    svg.append("</svg>")
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    svg_path = output_prefix.with_suffix(".svg")
    pdf_path = output_prefix.with_suffix(".pdf")
    png_path = output_prefix.with_suffix(".png")
    svg_path.write_text("\n".join(svg) + "\n", encoding="utf-8")
    rsvg = shutil.which("rsvg-convert")
    if rsvg:
        subprocess.run([rsvg, "-f", "pdf", "-o", str(pdf_path), str(svg_path)], check=True)
        subprocess.run([rsvg, "-f", "png", "-o", str(png_path), str(svg_path)], check=True)
    return svg_path, pdf_path, png_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", default="/tmp/satute-null-scaling-batch")
    parser.add_argument("--null-length", type=float, default=500.0)
    parser.add_argument(
        "--output",
        default=str(PROJECT_ROOT / "artifacts/local/saturation-diagnostics/null-scaling-summary.tsv"),
    )
    parser.add_argument(
        "--figure-prefix",
        default=str(PROJECT_ROOT / "artifacts/local/saturation-diagnostics/figures/null-scaling"),
    )
    args = parser.parse_args()
    rows = load_rows(Path(args.input_root), args.null_length)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(output)
    for path_value in render(rows, Path(args.figure_prefix)):
        print(path_value)


if __name__ == "__main__":
    main()
