#!/usr/bin/env python3

import argparse
import csv
import math
import subprocess
from pathlib import Path
from xml.sax.saxutils import escape


ROOT = Path(__file__).resolve().parents[2]
OUTDIR = ROOT / "artifacts" / "releases" / "generated-figures"

DATASETS = [
    {
        "label": "JC",
        "path": ROOT / "artifacts" / "cluster" / "lisc" / "paper-fig2-jc-reps1000" / "merged" / "head_to_head_summary.tsv",
        "simulation_model": "JC",
        "evaluation_model": "JC",
        "site_lengths": ["100", "1000", "10000"],
        "branch_ticks": [0.1, 0.5, 1.0, 2.5, 5.0, 10.0],
    },
    {
        "label": "GTR+F+G4",
        "path": ROOT
        / "artifacts"
        / "cluster"
        / "lisc"
        / "manuscript-ratehet-boundary-all-controls-reps1000"
        / "merged"
        / "head_to_head_summary.tsv",
        "simulation_model": "GTR_TOL_16S_G4",
        "evaluation_model": "GTR_TOL_16S_G4",
        "site_lengths": ["100", "250", "1000"],
        "branch_ticks": [4, 6, 8, 10, 12],
    },
    {
        "label": "LG+G4",
        "path": ROOT
        / "artifacts"
        / "cluster"
        / "lisc"
        / "manuscript-ratehet-boundary-all-controls-reps1000"
        / "merged"
        / "head_to_head_summary.tsv",
        "simulation_model": "LG_TOL_G4",
        "evaluation_model": "LG_TOL_G4",
        "site_lengths": ["100", "250", "1000"],
        "branch_ticks": [4, 6, 8, 10, 12],
    },
]

TREE_CASES = [
    ("five_external", "five-taxon external branch"),
    ("sixteen_internal", "16-taxon internal branch"),
]

FORMULAS = [
    ("dominant", "published SatuTe", "#0A3FBA"),
    ("eigenvalue_weighted", "eigenvalue-weighted", "#E69504"),
]

SITE_DASH = {
    "100": "8,5",
    "250": "4,4",
    "1000": "",
    "10000": "1,5",
}


def read_rows(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sx(value, xmin, xmax, left, width):
    return left + (math.log10(value) - math.log10(xmin)) / (math.log10(xmax) - math.log10(xmin)) * width


def sy(value, top, height):
    return top + (1.0 - value) * height


def path_data(points):
    if not points:
        return ""
    parts = [f"M {points[0][0]:.2f} {points[0][1]:.2f}"]
    parts.extend(f"L {x:.2f} {y:.2f}" for x, y in points[1:])
    return " ".join(parts)


def grouped(rows, dataset, tree_case):
    values = {}
    for row in rows:
        if row["simulation_model"] != dataset["simulation_model"]:
            continue
        if row["evaluation_model"] != dataset["evaluation_model"]:
            continue
        if row["scenario"] != "true_tree_fixed_lengths":
            continue
        if row["tree_case"] != tree_case:
            continue
        if row["fraction_informative"] == "":
            continue
        key = (row["formula"], row["nsites"])
        values.setdefault(key, {})[float(row["branch_length"])] = float(row["fraction_informative"])
    return values


def draw_panel(svg, rows, dataset, tree_case, title, left, top, width, height, show_y, show_x):
    data = grouped(rows, dataset, tree_case)
    all_branch_lengths = sorted({branch for series in data.values() for branch in series})
    xmin, xmax = min(all_branch_lengths), max(all_branch_lengths)

    svg.append(f'<text x="{left + width / 2:.1f}" y="{top - 18}" text-anchor="middle" font-size="17" font-weight="700">{escape(title)}</text>')
    svg.append(f'<rect x="{left}" y="{top}" width="{width}" height="{height}" fill="white" stroke="#444" stroke-width="1.1"/>')

    for value in [0, 0.25, 0.5, 0.75, 1.0]:
        y = sy(value, top, height)
        svg.append(f'<line x1="{left}" x2="{left + width}" y1="{y:.2f}" y2="{y:.2f}" stroke="#e6e6e6" stroke-width="1"/>')
        if show_y:
            svg.append(f'<text x="{left - 10}" y="{y + 5:.2f}" text-anchor="end" font-size="14">{value:g}</text>')

    for tick in dataset["branch_ticks"]:
        x = sx(tick, xmin, xmax, left, width)
        svg.append(f'<line x1="{x:.2f}" x2="{x:.2f}" y1="{top}" y2="{top + height}" stroke="#eeeeee" stroke-width="1"/>')
        if show_x:
            svg.append(f'<text x="{x:.2f}" y="{top + height + 22}" text-anchor="middle" font-size="14">{tick:g}</text>')

    for formula, _label, color in FORMULAS:
        for nsites in dataset["site_lengths"]:
            series = data.get((formula, nsites), {})
            points = [(sx(branch, xmin, xmax, left, width), sy(series[branch], top, height)) for branch in all_branch_lengths if branch in series]
            if not points:
                continue
            dash = SITE_DASH.get(nsites, "")
            dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
            svg.append(f'<path d="{path_data(points)}" fill="none" stroke="{color}" stroke-width="2.6"{dash_attr}/>')
            for x, y in points:
                svg.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3.1" fill="{color}" stroke="white" stroke-width="1"/>')

    if show_y:
        x = left - 58
        y = top + height / 2
        svg.append(f'<text x="{x}" y="{y:.1f}" text-anchor="middle" font-size="15" transform="rotate(-90 {x} {y:.1f})">fraction informative</text>')
    if show_x:
        svg.append(f'<text x="{left + width / 2:.1f}" y="{top + height + 50}" text-anchor="middle" font-size="15">target branch length</text>')


def legend(svg, x, y):
    svg.append(f'<text x="{x}" y="{y}" font-size="17" font-weight="700">Formula</text>')
    for i, (_formula, label, color) in enumerate(FORMULAS):
        yy = y + 28 + i * 28
        svg.append(f'<line x1="{x}" x2="{x + 42}" y1="{yy}" y2="{yy}" stroke="{color}" stroke-width="3.2"/>')
        svg.append(f'<text x="{x + 54}" y="{yy + 5}" font-size="15">{escape(label)}</text>')

    y2 = y + 100
    svg.append(f'<text x="{x}" y="{y2}" font-size="17" font-weight="700">Alignment length</text>')
    for i, nsites in enumerate(["100", "250", "1000", "10000"]):
        yy = y2 + 28 + i * 28
        dash = SITE_DASH.get(nsites, "")
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        svg.append(f'<line x1="{x}" x2="{x + 42}" y1="{yy}" y2="{yy}" stroke="#222222" stroke-width="2.8"{dash_attr}/>')
        svg.append(f'<text x="{x + 54}" y="{yy + 5}" font-size="15">{nsites} sites</text>')


def main():
    parser = argparse.ArgumentParser(description="Render the reviewed relative-weighting power summary.")
    parser.add_argument("--output-dir", type=Path, default=OUTDIR)
    parser.add_argument("--jc-summary", type=Path, default=DATASETS[0]["path"])
    parser.add_argument("--rate-summary", type=Path, default=DATASETS[1]["path"])
    args = parser.parse_args()
    datasets = [dict(dataset) for dataset in DATASETS]
    datasets[0]["path"] = args.jc_summary
    datasets[1]["path"] = args.rate_summary
    datasets[2]["path"] = args.rate_summary
    out_svg = args.output_dir / "figure_power_curve_summary.svg"
    out_pdf = args.output_dir / "figure_power_curve_summary.pdf"

    rows_by_path = {}
    for dataset in datasets:
        rows_by_path.setdefault(dataset["path"], read_rows(dataset["path"]))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    width, height = 1480, 1300
    panel_w, panel_h = 395, 220
    x0, y0 = 135, 130
    xgap, ygap = 105, 170

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{width / 2:.1f}" y="46" text-anchor="middle" font-size="30" font-weight="700">Full-alignment simulation curves</text>',
        f'<text x="{width / 2:.1f}" y="78" text-anchor="middle" font-size="17" fill="#333333">True-tree fixed-length scenario; 1,000 replicates per point</text>',
    ]

    for row_idx, dataset in enumerate(datasets):
        top = y0 + row_idx * (panel_h + ygap)
        svg.append(f'<text x="30" y="{top + panel_h / 2:.1f}" text-anchor="middle" font-size="21" font-weight="700" transform="rotate(-90 30 {top + panel_h / 2:.1f})">{escape(dataset["label"])}</text>')
        rows = rows_by_path[dataset["path"]]
        for col_idx, (tree_case, tree_label) in enumerate(TREE_CASES):
            left = x0 + col_idx * (panel_w + xgap)
            draw_panel(
                svg,
                rows,
                dataset,
                tree_case,
                tree_label,
                left,
                top,
                panel_w,
                panel_h,
                show_y=(col_idx == 0),
                show_x=True,
            )

    legend(svg, 1160, 140)
    svg.append('<text x="135" y="1250" font-size="14" fill="#555555">JC panels use 100, 1,000 and 10,000 sites; GTR+F+G4 and LG+G4 panels use 100, 250 and 1,000 sites.</text>')
    svg.append("</svg>")

    out_svg.write_text("\n".join(svg) + "\n", encoding="utf-8")
    subprocess.run(["rsvg-convert", "-f", "pdf", "-o", str(out_pdf), str(out_svg)], check=True)
    print(out_svg)
    print(out_pdf)


if __name__ == "__main__":
    main()
