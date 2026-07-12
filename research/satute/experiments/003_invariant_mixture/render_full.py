#!/usr/bin/env python3
"""Render the dense full +I analysis as publication-quality vector panels."""

import argparse
import csv
import math
import shutil
import subprocess
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]


METHOD_STYLE = {
    "legacy_counterfactual": ("Legacy counterfactual", "#D97706", "8,5"),
    "refined_invariant_null": ("Refined invariant null", "#0A3FBA", ""),
}
MODEL_DASH = {"+I+G4": "", "+I+R4": "3,3"}
SITE_LENGTHS = (100, 500, 1000, 5000, 10000)


def read_rows(path, text_fields):
    with path.open("r", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        for key, value in list(row.items()):
            if key not in text_fields:
                row[key] = float(value)
    return rows


def log_scale(value, lower, upper, start, length):
    return start + length * (math.log10(value) - math.log10(lower)) / (
        math.log10(upper) - math.log10(lower)
    )


def y_scale(value, lower, upper, start, length):
    return start + length * (1.0 - (value - lower) / (upper - lower))


def line_axes(svg, left, top, width, height, title, subtitle, xmin, xmax, ymin, ymax, ylabel):
    svg.append(
        f'<text x="{left}" y="{top-31}" font-size="16" font-weight="700" '
        f'fill="#202124">{title}</text>'
    )
    svg.append(
        f'<text x="{left}" y="{top-11}" font-size="10.5" fill="#5f6368">{subtitle}</text>'
    )
    svg.append(
        f'<rect x="{left}" y="{top}" width="{width}" height="{height}" '
        'fill="#fff" stroke="#3c4043"/>'
    )
    for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
        value = ymin + fraction * (ymax - ymin)
        y = y_scale(value, ymin, ymax, top, height)
        svg.append(
            f'<line x1="{left}" x2="{left+width}" y1="{y:.2f}" y2="{y:.2f}" '
            'stroke="#eceff1"/>'
        )
        svg.append(
            f'<text x="{left-8}" y="{y+4:.2f}" text-anchor="end" font-size="9" '
            f'fill="#5f6368">{value:.3g}</text>'
        )
    ticks = (0.5, 1, 2, 5, 10, 20, 50, 100, 200, 1000, 10000)
    for tick in ticks:
        if xmin <= tick <= xmax:
            x = log_scale(tick, xmin, xmax, left, width)
            svg.append(
                f'<text x="{x:.2f}" y="{top+height+17}" text-anchor="middle" '
                f'font-size="8" fill="#5f6368">{tick:g}</text>'
            )
    svg.append(
        f'<text x="{left+width/2}" y="{top+height+38}" text-anchor="middle" '
        'font-size="10" fill="#3c4043">target branch length, t (log scale)</text>'
    )
    svg.append(
        f'<text x="{left-56}" y="{top+height/2}" text-anchor="middle" '
        f'font-size="10" fill="#3c4043" transform="rotate(-90 {left-56} '
        f'{top+height/2})">{ylabel}</text>'
    )


def draw_line_band(svg, rows, left, top, width, height, xmin, xmax, ymin, ymax, center, low, high, color, dash):
    rows = sorted(rows, key=lambda row: row["branch_length"])
    upper = [
        (
            log_scale(row["branch_length"], xmin, xmax, left, width),
            y_scale(row[high], ymin, ymax, top, height),
        )
        for row in rows
    ]
    lower = [
        (
            log_scale(row["branch_length"], xmin, xmax, left, width),
            y_scale(row[low], ymin, ymax, top, height),
        )
        for row in reversed(rows)
    ]
    polygon = " ".join(f"{x:.2f},{y:.2f}" for x, y in upper + lower)
    svg.append(f'<polygon points="{polygon}" fill="{color}" fill-opacity="0.10"/>')
    points = [
        (
            log_scale(row["branch_length"], xmin, xmax, left, width),
            y_scale(row[center], ymin, ymax, top, height),
        )
        for row in rows
    ]
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    svg.append(
        '<polyline points="'
        + " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
        + f'" fill="none" stroke="{color}" stroke-width="2.4"{dash_attr}/>'
    )
    for x, y in points:
        svg.append(
            f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2.4" fill="#fff" '
            f'stroke="{color}" stroke-width="1.4"/>'
        )


def interpolate(start, end, fraction):
    fraction = max(0.0, min(1.0, fraction))
    rgb = tuple(round(a + fraction * (b - a)) for a, b in zip(start, end))
    return "#" + "".join(f"{value:02x}" for value in rgb)


def signed_color(value, bound):
    white = (255, 255, 255)
    if value >= 0.0:
        return interpolate(white, (10, 63, 186), value / bound)
    return interpolate(white, (182, 73, 125), -value / bound)


def branch_boundaries(branches):
    logs = [math.log10(value) for value in branches]
    boundaries = [logs[0] - (logs[1] - logs[0]) / 2.0]
    boundaries.extend((left + right) / 2.0 for left, right in zip(logs[:-1], logs[1:]))
    boundaries.append(logs[-1] + (logs[-1] - logs[-2]) / 2.0)
    return [10**value for value in boundaries]


def heatmap(svg, rows, left, top, width, height, title, subtitle, bound):
    branches = sorted({row["branch_length"] for row in rows})
    boundaries = branch_boundaries(branches)
    xmin, xmax = boundaries[0], boundaries[-1]
    svg.append(
        f'<text x="{left}" y="{top-31}" font-size="16" font-weight="700" '
        f'fill="#202124">{title}</text>'
    )
    svg.append(
        f'<text x="{left}" y="{top-11}" font-size="10.5" fill="#5f6368">{subtitle}</text>'
    )
    cell_h = height / len(SITE_LENGTHS)
    lookup = {(row["nsites"], row["branch_length"]): row for row in rows}
    for row_index, nsites in enumerate(SITE_LENGTHS):
        y = top + row_index * cell_h
        for branch_index, branch in enumerate(branches):
            x0 = log_scale(boundaries[branch_index], xmin, xmax, left, width)
            x1 = log_scale(boundaries[branch_index + 1], xmin, xmax, left, width)
            value = lookup[(nsites, branch)]["refined_minus_legacy"]
            svg.append(
                f'<rect x="{x0:.2f}" y="{y:.2f}" width="{x1-x0:.2f}" '
                f'height="{cell_h:.2f}" fill="{signed_color(value, bound)}" '
                'stroke="#ffffff" stroke-width="0.45"/>'
            )
        svg.append(
            f'<text x="{left-8}" y="{y+cell_h/2+4:.2f}" text-anchor="end" '
            f'font-size="9" fill="#5f6368">{nsites:,}</text>'
        )
    svg.append(
        f'<rect x="{left}" y="{top}" width="{width}" height="{height}" '
        'fill="none" stroke="#3c4043"/>'
    )
    for tick in (0.5, 1, 2, 5, 10, 20, 50, 100, 200, 1000, 10000):
        if branches[0] <= tick <= branches[-1]:
            x = log_scale(tick, xmin, xmax, left, width)
            svg.append(
                f'<text x="{x:.2f}" y="{top+height+17}" text-anchor="middle" '
                f'font-size="8" fill="#5f6368">{tick:g}</text>'
            )
    svg.append(
        f'<text x="{left+width/2}" y="{top+height+38}" text-anchor="middle" '
        'font-size="10" fill="#3c4043">target branch length, t (log scale)</text>'
    )
    svg.append(
        f'<text x="{left-56}" y="{top+height/2}" text-anchor="middle" '
        f'font-size="10" fill="#3c4043" transform="rotate(-90 {left-56} '
        f'{top+height/2})">alignment length, n</text>'
    )
    # Diverging color scale.
    bar_x, bar_y, bar_w, bar_h = left + width - 150, top - 29, 145, 8
    segments = 80
    for index in range(segments):
        value = -bound + (2.0 * bound) * index / (segments - 1)
        svg.append(
            f'<rect x="{bar_x+bar_w*index/segments:.2f}" y="{bar_y}" '
            f'width="{bar_w/segments+0.3:.2f}" height="{bar_h}" '
            f'fill="{signed_color(value, bound)}"/>'
        )
    svg.append(
        f'<text x="{bar_x}" y="{bar_y-3}" font-size="7.5" fill="#5f6368">{-bound:.3f}</text>'
    )
    svg.append(
        f'<text x="{bar_x+bar_w/2}" y="{bar_y-3}" text-anchor="middle" '
        'font-size="7.5" fill="#5f6368">0</text>'
    )
    svg.append(
        f'<text x="{bar_x+bar_w}" y="{bar_y-3}" text-anchor="end" '
        f'font-size="7.5" fill="#5f6368">+{bound:.3f}</text>'
    )


def n_axes(svg, left, top, width, height, title, subtitle, ymin, ymax, ylabel):
    svg.append(
        f'<text x="{left}" y="{top-31}" font-size="16" font-weight="700" '
        f'fill="#202124">{title}</text>'
    )
    svg.append(
        f'<text x="{left}" y="{top-11}" font-size="10.5" fill="#5f6368">{subtitle}</text>'
    )
    svg.append(
        f'<rect x="{left}" y="{top}" width="{width}" height="{height}" '
        'fill="#fff" stroke="#3c4043"/>'
    )
    for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
        value = ymin + fraction * (ymax - ymin)
        y = y_scale(value, ymin, ymax, top, height)
        svg.append(
            f'<line x1="{left}" x2="{left+width}" y1="{y:.2f}" y2="{y:.2f}" '
            'stroke="#eceff1"/>'
        )
        svg.append(
            f'<text x="{left-8}" y="{y+4:.2f}" text-anchor="end" font-size="9" '
            f'fill="#5f6368">{value:.3g}</text>'
        )
    for tick in SITE_LENGTHS:
        x = log_scale(tick, 100, 10000, left, width)
        svg.append(
            f'<text x="{x:.2f}" y="{top+height+17}" text-anchor="middle" '
            f'font-size="8" fill="#5f6368">{tick:g}</text>'
        )
    svg.append(
        f'<text x="{left+width/2}" y="{top+height+38}" text-anchor="middle" '
        'font-size="10" fill="#3c4043">alignment length, n (log scale)</text>'
    )
    svg.append(
        f'<text x="{left-56}" y="{top+height/2}" text-anchor="middle" '
        f'font-size="10" fill="#3c4043" transform="rotate(-90 {left-56} '
        f'{top+height/2})">{ylabel}</text>'
    )


def simple_n_line(svg, rows, left, top, width, height, ymin, ymax, value_key, color, dash):
    rows = sorted(rows, key=lambda row: row["nsites"])
    points = [
        (
            log_scale(row["nsites"], 100, 10000, left, width),
            y_scale(row[value_key], ymin, ymax, top, height),
        )
        for row in rows
    ]
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    svg.append(
        '<polyline points="'
        + " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
        + f'" fill="none" stroke="{color}" stroke-width="2.3"{dash_attr}/>'
    )
    for x, y in points:
        svg.append(
            f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2.6" fill="#fff" '
            f'stroke="{color}" stroke-width="1.4"/>'
        )


def render(summary, paired, output_prefix, png_scale):
    width, height = 1600, 1400
    panel_w, panel_h = 620, 285
    lefts = (115, 865)
    tops = (175, 590, 1005)
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="800" y="40" text-anchor="middle" font-size="27" font-weight="700" fill="#202124">Full paired invariant-mixture analysis</text>',
        '<text x="800" y="68" text-anchor="middle" font-size="12.5" fill="#5f6368">Dense branch grid; 20,000 null and 20,000 alternative replicates per model × branch × alignment-length cell</text>',
    ]
    branches = sorted({row["branch_length"] for row in summary})
    xmin, xmax = min(branches), max(branches)

    for index, model in enumerate(("+I+G4", "+I+R4")):
        left, top = lefts[index], tops[0]
        line_axes(
            svg,
            left,
            top,
            panel_w,
            panel_h,
            f"{chr(65+index)}  {model} calibrated power",
            "n = 1,000; bands are Wilson 95% intervals",
            xmin,
            xmax,
            0.0,
            1.0,
            "probability informative",
        )
        for method, (_label, color, dash) in METHOD_STYLE.items():
            rows = [
                row
                for row in summary
                if row["model"] == model
                and row["method"] == method
                and row["nsites"] == 1000
            ]
            draw_line_band(
                svg,
                rows,
                left,
                top,
                panel_w,
                panel_h,
                xmin,
                xmax,
                0.0,
                1.0,
                "calibrated_power",
                "calibrated_power_ci_low",
                "calibrated_power_ci_high",
                color,
                dash,
            )

    heat_rows = [
        row
        for row in paired
        if row["sample"] == "alternative" and row["rule"] == "calibrated"
    ]
    bound = max(abs(row["refined_minus_legacy"]) for row in heat_rows)
    bound = math.ceil(bound * 1000.0) / 1000.0
    for index, model in enumerate(("+I+G4", "+I+R4")):
        rows = [row for row in heat_rows if row["model"] == model]
        heatmap(
            svg,
            rows,
            lefts[index],
            tops[1],
            panel_w,
            panel_h,
            f"{chr(67+index)}  {model} paired power difference",
            "refined minus legacy after independent null calibration",
            bound,
        )

    # Extreme-branch null behavior.
    left, top = lefts[0], tops[2]
    n_axes(
        svg,
        left,
        top,
        panel_w,
        panel_h,
        "E  Extreme-branch nominal false positives",
        "t = 10,000 under the correct invariant-mixture null",
        0.0,
        1.0,
        "rejection fraction",
    )
    y05 = y_scale(0.05, 0.0, 1.0, top, panel_h)
    svg.append(
        f'<line x1="{left}" x2="{left+panel_w}" y1="{y05:.2f}" y2="{y05:.2f}" '
        'stroke="#3c4043" stroke-dasharray="4,4"/>'
    )
    for model in ("+I+G4", "+I+R4"):
        for method, (_label, color, method_dash) in METHOD_STYLE.items():
            rows = [
                row
                for row in summary
                if row["model"] == model
                and row["method"] == method
                and row["branch_length"] == 10000
            ]
            dash = MODEL_DASH[model] or method_dash
            if MODEL_DASH[model] and method_dash:
                dash = "9,3,2,3"
            simple_n_line(
                svg,
                rows,
                left,
                top,
                panel_w,
                panel_h,
                0.0,
                1.0,
                "null_nominal_rejection",
                color,
                dash,
            )

    left, top = lefts[1], tops[2]
    extreme = [row for row in summary if row["branch_length"] == 10000]
    ymax = max(row["empirical_threshold"] for row in extreme) * 1.05
    n_axes(
        svg,
        left,
        top,
        panel_w,
        panel_h,
        "F  Extreme-branch empirical thresholds",
        "95th percentile trained on an independent half-null sample",
        0.0,
        ymax,
        "calibrated Z threshold",
    )
    for model in ("+I+G4", "+I+R4"):
        for method, (_label, color, method_dash) in METHOD_STYLE.items():
            rows = [
                row
                for row in extreme
                if row["model"] == model and row["method"] == method
            ]
            dash = MODEL_DASH[model] or method_dash
            if MODEL_DASH[model] and method_dash:
                dash = "9,3,2,3"
            simple_n_line(
                svg,
                rows,
                left,
                top,
                panel_w,
                panel_h,
                0.0,
                ymax,
                "empirical_threshold",
                color,
                dash,
            )

    # Legends: color identifies formula, line texture identifies rate mixture.
    legend_y = 108
    for index, (_method, (label, color, dash)) in enumerate(METHOD_STYLE.items()):
        x = 280 + index * 350
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        svg.append(
            f'<line x1="{x}" x2="{x+44}" y1="{legend_y}" y2="{legend_y}" '
            f'stroke="{color}" stroke-width="2.7"{dash_attr}/>'
        )
        svg.append(
            f'<text x="{x+54}" y="{legend_y+4}" font-size="11" fill="#3c4043">{label}</text>'
        )
    x = 1030
    svg.append(
        f'<line x1="{x}" x2="{x+44}" y1="{legend_y}" y2="{legend_y}" '
        'stroke="#3c4043" stroke-width="2.3"/>'
    )
    svg.append(
        f'<text x="{x+52}" y="{legend_y+4}" font-size="11" fill="#3c4043">+I+G4</text>'
    )
    x = 1200
    svg.append(
        f'<line x1="{x}" x2="{x+44}" y1="{legend_y}" y2="{legend_y}" '
        'stroke="#3c4043" stroke-width="2.3" stroke-dasharray="3,3"/>'
    )
    svg.append(
        f'<text x="{x+52}" y="{legend_y+4}" font-size="11" fill="#3c4043">+I+R4</text>'
    )
    svg.append(
        '<text x="115" y="1380" font-size="10.5" fill="#5f6368">Vector PDF/SVG and 3× PNG exports. Heatmaps use one symmetric scale across both models; blue favors the refined statistic and pink favors the legacy statistic.</text>'
    )
    svg.append("</svg>")

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    svg_path = output_prefix.with_suffix(".svg")
    svg_path.write_text("\n".join(svg) + "\n", encoding="utf-8")
    rsvg = shutil.which("rsvg-convert")
    if rsvg:
        subprocess.run(
            [
                rsvg,
                "-f",
                "png",
                "-z",
                str(png_scale),
                "-o",
                str(output_prefix.with_suffix(".png")),
                str(svg_path),
            ],
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
        default=str(PROJECT_ROOT / "artifacts/local/weighted-invariant-full"),
    )
    parser.add_argument(
        "--output-prefix",
        default=str(PROJECT_ROOT / "artifacts/local/weighted-invariant-full/figures/full-analysis"),
    )
    parser.add_argument("--png-scale", type=float, default=3.0)
    args = parser.parse_args()
    root = Path(args.input_dir)
    summary = read_rows(
        root / "weighted_invariant_summary.tsv", {"model", "method"}
    )
    paired = read_rows(
        root / "weighted_invariant_paired.tsv", {"model", "rule", "sample"}
    )
    render(summary, paired, Path(args.output_prefix), args.png_scale)
    print(Path(args.output_prefix).with_suffix(".svg"))


if __name__ == "__main__":
    main()
