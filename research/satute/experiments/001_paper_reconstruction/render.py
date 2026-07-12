#!/usr/bin/env python3

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


SCENARIOS = [
    ("true_tree_fixed_lengths", "#d95f8d", "true tree, fixed lengths"),
    ("true_topology_ml_lengths", "#7b3294", "true topology, ML lengths"),
    ("ml_tree_unadjusted", "#2c7fb8", "ML tree, unadjusted"),
    ("ml_tree_bonferroni", "#1a9850", "ML tree, Bonferroni"),
]

SITE_STYLES = {
    "100": "6,4",
    "1000": "",
    "10000": "1,4",
}


def load_rows(path):
    with open(path, "r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sx(value, xmin, xmax, left, width):
    return left + (math.log10(value) - math.log10(xmin)) / (math.log10(xmax) - math.log10(xmin)) * width


def sy(value, top, height):
    return top + (1.0 - value) * height


def line_path(points):
    if not points:
        return ""
    parts = [f"M {points[0][0]:.2f} {points[0][1]:.2f}"]
    for x, y in points[1:]:
        parts.append(f"L {x:.2f} {y:.2f}")
    return " ".join(parts)


def panel(svg, rows, tree_case, title, left, top, width, height):
    subset = [row for row in rows if row["tree_case"] == tree_case and row["evaluated"] != "0"]
    branch_lengths = sorted({float(row["branch_length"]) for row in subset})
    xmin, xmax = min(branch_lengths), max(branch_lengths)

    svg.append(f'<g transform="translate(0,0)">')
    svg.append(f'<text x="{left + width / 2:.1f}" y="{top - 18}" text-anchor="middle" font-size="16" font-weight="600">{title}</text>')
    svg.append(f'<rect x="{left}" y="{top}" width="{width}" height="{height}" fill="white" stroke="#333" stroke-width="1"/>')

    for ytick in [0.0, 0.25, 0.5, 0.75, 1.0]:
        y = sy(ytick, top, height)
        svg.append(f'<line x1="{left}" x2="{left + width}" y1="{y:.2f}" y2="{y:.2f}" stroke="#e8e8e8" stroke-width="1"/>')
        svg.append(f'<text x="{left - 8}" y="{y + 4:.2f}" text-anchor="end" font-size="11">{ytick:.2g}</text>')

    for xtick in branch_lengths:
        x = sx(xtick, xmin, xmax, left, width)
        svg.append(f'<line x1="{x:.2f}" x2="{x:.2f}" y1="{top}" y2="{top + height}" stroke="#f0f0f0" stroke-width="1"/>')
        svg.append(f'<text x="{x:.2f}" y="{top + height + 18}" text-anchor="middle" font-size="11">{xtick:g}</text>')

    grouped = defaultdict(dict)
    for row in subset:
        grouped[(row["scenario"], row["nsites"])][float(row["branch_length"])] = float(row["fraction_informative"])

    for scenario, color, _label in SCENARIOS:
        for nsites in sorted({row["nsites"] for row in subset}, key=int):
            values = grouped.get((scenario, nsites), {})
            points = [
                (sx(blen, xmin, xmax, left, width), sy(values[blen], top, height))
                for blen in branch_lengths
                if blen in values
            ]
            dash = SITE_STYLES.get(nsites, "")
            dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
            svg.append(
                f'<path d="{line_path(points)}" fill="none" stroke="{color}" stroke-width="2.4"{dash_attr}/>'
            )
            for x, y in points:
                svg.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3" fill="{color}" stroke="white" stroke-width="1"/>')

    svg.append(f'<text x="{left + width / 2:.1f}" y="{top + height + 42}" text-anchor="middle" font-size="12">branch length AB, log scale</text>')
    svg.append(f'<text x="{left - 48}" y="{top + height / 2:.1f}" text-anchor="middle" font-size="12" transform="rotate(-90 {left - 48} {top + height / 2:.1f})">fraction informative</text>')
    svg.append("</g>")


def legend(svg, x, y):
    svg.append(f'<text x="{x}" y="{y}" font-size="14" font-weight="600">Scenario</text>')
    dy = 22
    for i, (_scenario, color, label) in enumerate(SCENARIOS):
        yy = y + 22 + i * dy
        svg.append(f'<line x1="{x}" x2="{x + 34}" y1="{yy}" y2="{yy}" stroke="{color}" stroke-width="3"/>')
        svg.append(f'<text x="{x + 44}" y="{yy + 4}" font-size="12">{label}</text>')
    y2 = y + 22 + len(SCENARIOS) * dy + 14
    svg.append(f'<text x="{x}" y="{y2}" font-size="14" font-weight="600">Sequence length</text>')
    for i, (nsites, dash) in enumerate([("100", "6,4"), ("1000", "")]):
        yy = y2 + 22 + i * dy
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        svg.append(f'<line x1="{x}" x2="{x + 34}" y1="{yy}" y2="{yy}" stroke="#222" stroke-width="2.4"{dash_attr}/>')
        svg.append(f'<text x="{x + 44}" y="{yy + 4}" font-size="12">{nsites} sites</text>')


def main():
    project = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description="Render the reconstructed SatuTe simulation curves.")
    parser.add_argument(
        "--summary",
        type=Path,
        default=project / "artifacts" / "releases" / "legacy" / "fig2_reduced_summary.tsv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=project / "artifacts" / "local" / "paper-reconstruction" / "fig2_reduced_reconstruction.svg",
    )
    args = parser.parse_args()
    summary = args.summary
    output = args.output
    rows = load_rows(summary)

    svg = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1160" height="470" viewBox="0 0 1160 470">',
        '<rect width="1160" height="470" fill="#ffffff"/>',
        '<text x="580" y="28" text-anchor="middle" font-size="19" font-weight="700">Reduced reconstruction of the SatuTe simulation design</text>',
    ]
    panel(svg, rows, "five_external", "a) five-taxon external branch", 78, 78, 390, 250)
    panel(svg, rows, "sixteen_internal", "b) 16-taxon internal branch", 555, 78, 390, 250)
    legend(svg, 975, 92)
    svg.append('<text x="78" y="445" font-size="11" fill="#555">Reduced run: 2 replicates per point; branch lengths 0.1, 1.0, 5.0; site lengths 100 and 1000. Use the reconstruction script with --paper-grid --reps 1000 for the full paper-scale grid.</text>')
    svg.append("</svg>")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(svg) + "\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
