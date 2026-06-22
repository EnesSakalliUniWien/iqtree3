#!/usr/bin/env python3

import argparse
import csv
import math
import subprocess
from collections import defaultdict
from pathlib import Path


PAPER_BRANCH_LENGTHS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.5, 10.0]
PAPER_SITE_LENGTHS = ["100", "1000", "10000"]
GTR_PF06346_MODEL = "GTR{0.6676,3.7807,4.2833,0.5354,0.8718,1.0}+F{0.125,0.436,0.191,0.245}"
MODEL_ALIASES = {
    "JC": "JC",
    "GTR_PF06346": GTR_PF06346_MODEL,
    "GTR_EvoNAPS_PF06346": GTR_PF06346_MODEL,
}
TREE_CASES_ALL = [
    ("five_external", "a) five-taxon external branch"),
    ("sixteen_internal", "b) 16-taxon internal branch"),
]
SCENARIOS = [
    ("true_tree_fixed_lengths", "#d95f8d", "true tree, fixed lengths"),
    ("true_topology_ml_lengths", "#7b3294", "true topology, ML lengths"),
    ("ml_tree_unadjusted", "#2c7fb8", "ML tree, unadjusted"),
    ("ml_tree_bonferroni", "#1a9850", "ML tree, Bonferroni"),
]
FORMULAS = [
    ("dominant", "published dominant coefficient"),
    ("all_unweighted", "unweighted all-mode coefficient"),
    ("eigenvalue_weighted", "eigenvalue-weighted coefficient"),
]
SITE_STYLES = {
    "100": "6,4",
    "1000": "",
    "10000": "1,4",
}


def load_rows(path):
    with open(path, "r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def float_key(value):
    return round(float(value), 10)


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


def selected_tree_cases(value):
    labels = dict(TREE_CASES_ALL)
    tree_cases = []
    for tree_case in [item.strip() for item in value.split(",") if item.strip()]:
        if tree_case not in labels:
            raise SystemExit(f"Unknown tree case {tree_case!r}. Expected one of: {', '.join(labels)}")
        tree_cases.append((tree_case, labels[tree_case]))
    if not tree_cases:
        raise SystemExit("--tree-cases must contain at least one tree case")
    return tree_cases


def resolve_model(value, rows):
    models = sorted({row["model"] for row in rows})
    if not value:
        if len(models) == 1:
            return models[0]
        raise SystemExit(f"Summary contains multiple models; pass --model. Models: {', '.join(models)}")
    value = MODEL_ALIASES.get(value, value)
    if value not in models:
        raise SystemExit(f"Model {value!r} is not present in the summary. Models: {', '.join(models)}")
    return value


def model_label(model):
    if model == GTR_PF06346_MODEL:
        return "GTR_PF06346"
    return "".join(ch if ch.isalnum() else "_" for ch in model).strip("_")


def validate_full(rows, expected_reps, allow_incomplete, tree_cases):
    problems = []
    by_key = {
        (
            row["tree_case"],
            row["formula"],
            row["scenario"],
            row["nsites"],
            float_key(row["branch_length"]),
        ): row
        for row in rows
    }
    expected_branch_lengths = {float_key(value) for value in PAPER_BRANCH_LENGTHS}
    for tree_case, _ in tree_cases:
        for formula, _ in FORMULAS:
            for scenario, _color, _label in SCENARIOS:
                for nsites in PAPER_SITE_LENGTHS:
                    seen = {
                        branch
                        for (tc, ff, sc, ns, branch), _row in by_key.items()
                        if tc == tree_case and ff == formula and sc == scenario and ns == nsites
                    }
                    missing_branch_lengths = sorted(expected_branch_lengths - seen)
                    if missing_branch_lengths:
                        problems.append(
                            f"{tree_case}/{formula}/{scenario}/n{nsites} missing branch lengths {missing_branch_lengths}"
                        )
                    for branch_length in expected_branch_lengths & seen:
                        row = by_key[(tree_case, formula, scenario, nsites, branch_length)]
                        evaluated = int(row["evaluated"])
                        missing = int(row["missing_split"])
                        total = evaluated + missing
                        if total != expected_reps:
                            problems.append(
                                f"{tree_case}/{formula}/{scenario}/n{nsites}/b{branch_length:g} has "
                                f"evaluated+missing={total}, expected {expected_reps}"
                            )
                        if scenario in {"true_tree_fixed_lengths", "true_topology_ml_lengths"} and missing != 0:
                            problems.append(
                                f"{tree_case}/{formula}/{scenario}/n{nsites}/b{branch_length:g} has missing splits"
                            )

    if problems and not allow_incomplete:
        sample = "\n".join(f"- {problem}" for problem in problems[:20])
        extra = "" if len(problems) <= 20 else f"\n... {len(problems) - 20} more"
        raise SystemExit(
            "Input is not a complete paper-scale head-to-head run. "
            "Use --allow-incomplete only for plotting development artifacts.\n"
            f"{sample}{extra}"
        )
    return problems


def grouped_values(rows, formula, tree_case):
    grouped = defaultdict(dict)
    for row in rows:
        if row["tree_case"] != tree_case or row["formula"] != formula:
            continue
        if row["evaluated"] == "0" or row["fraction_informative"] == "":
            continue
        grouped[(row["scenario"], row["nsites"])][float(row["branch_length"])] = float(row["fraction_informative"])
    return grouped


def panel(svg, rows, formula, tree_case, title, left, top, width, height, show_y_label=True, show_x_label=True):
    grouped = grouped_values(rows, formula, tree_case)
    branch_lengths = PAPER_BRANCH_LENGTHS
    xmin, xmax = min(branch_lengths), max(branch_lengths)

    svg.append(f'<text x="{left + width / 2:.1f}" y="{top - 16}" text-anchor="middle" font-size="15" font-weight="600">{title}</text>')
    svg.append(f'<rect x="{left}" y="{top}" width="{width}" height="{height}" fill="white" stroke="#333" stroke-width="1"/>')

    for ytick in [0.0, 0.25, 0.5, 0.75, 1.0]:
        y = sy(ytick, top, height)
        svg.append(f'<line x1="{left}" x2="{left + width}" y1="{y:.2f}" y2="{y:.2f}" stroke="#e8e8e8" stroke-width="1"/>')
        if show_y_label:
            svg.append(f'<text x="{left - 8}" y="{y + 4:.2f}" text-anchor="end" font-size="10">{ytick:.2g}</text>')

    for xtick in branch_lengths:
        x = sx(xtick, xmin, xmax, left, width)
        svg.append(f'<line x1="{x:.2f}" x2="{x:.2f}" y1="{top}" y2="{top + height}" stroke="#f0f0f0" stroke-width="1"/>')
        if xtick in {0.1, 0.5, 1.0, 2.5, 5.0, 10.0}:
            svg.append(f'<text x="{x:.2f}" y="{top + height + 16}" text-anchor="middle" font-size="10">{xtick:g}</text>')

    for scenario, color, _label in SCENARIOS:
        for nsites in PAPER_SITE_LENGTHS:
            values = grouped.get((scenario, nsites), {})
            points = [
                (sx(blen, xmin, xmax, left, width), sy(values[blen], top, height))
                for blen in branch_lengths
                if blen in values
            ]
            if not points:
                continue
            dash = SITE_STYLES.get(nsites, "")
            dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
            svg.append(f'<path d="{line_path(points)}" fill="none" stroke="{color}" stroke-width="2.1"{dash_attr}/>')
            for x, y in points:
                svg.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2.5" fill="{color}" stroke="white" stroke-width="0.8"/>')

    if show_x_label:
        svg.append(f'<text x="{left + width / 2:.1f}" y="{top + height + 38}" text-anchor="middle" font-size="11">branch length AB, log scale</text>')
    if show_y_label:
        svg.append(f'<text x="{left - 42}" y="{top + height / 2:.1f}" text-anchor="middle" font-size="11" transform="rotate(-90 {left - 42} {top + height / 2:.1f})">fraction informative</text>')


def legend(svg, x, y):
    svg.append(f'<text x="{x}" y="{y}" font-size="14" font-weight="600">Scenario</text>')
    dy = 21
    for i, (_scenario, color, label) in enumerate(SCENARIOS):
        yy = y + 21 + i * dy
        svg.append(f'<line x1="{x}" x2="{x + 32}" y1="{yy}" y2="{yy}" stroke="{color}" stroke-width="3"/>')
        svg.append(f'<text x="{x + 42}" y="{yy + 4}" font-size="11">{label}</text>')
    y2 = y + 21 + len(SCENARIOS) * dy + 12
    svg.append(f'<text x="{x}" y="{y2}" font-size="14" font-weight="600">Sequence length</text>')
    for i, nsites in enumerate(PAPER_SITE_LENGTHS):
        yy = y2 + 21 + i * dy
        dash = SITE_STYLES.get(nsites, "")
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        svg.append(f'<line x1="{x}" x2="{x + 32}" y1="{yy}" y2="{yy}" stroke="#222" stroke-width="2.4"{dash_attr}/>')
        svg.append(f'<text x="{x + 42}" y="{yy + 4}" font-size="11">{nsites} sites</text>')


def write_baseline_figure(rows, output, problems, tree_cases, model):
    panel_count = len(tree_cases)
    width = 720 if panel_count == 1 else 1160
    panel_width = 430 if panel_count == 1 else 390
    panel_lefts = [92] if panel_count == 1 else [78, 555]
    legend_x = 555 if panel_count == 1 else 975
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="470" viewBox="0 0 {width} 470">',
        f'<rect width="{width}" height="470" fill="#ffffff"/>',
        f'<text x="{width / 2:.1f}" y="28" text-anchor="middle" font-size="19" font-weight="700">SatuTe simulation reproduction: published dominant statistic ({model_label(model)})</text>',
    ]
    for index, (tree_case, tree_label) in enumerate(tree_cases):
        panel(svg, rows, "dominant", tree_case, tree_label, panel_lefts[index], 78, panel_width, 250)
    legend(svg, legend_x, 92)
    note = "Full paper grid: 1,000 replicates per point; branch lengths 0.1 to 10; site lengths 100, 1,000 and 10,000."
    if problems:
        note = "Incomplete plotting-development output; not for manuscript use."
    svg.append(f'<text x="78" y="445" font-size="11" fill="#555">{note}</text>')
    svg.append("</svg>")
    output.write_text("\n".join(svg) + "\n", encoding="utf-8")


def write_formula_grid(rows, output, problems, tree_cases, model):
    width = 1420
    height = 540 if len(tree_cases) == 1 else 880
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<rect width="{width}" height="{height}" fill="#ffffff"/>',
        f'<text x="710" y="28" text-anchor="middle" font-size="19" font-weight="700">Head-to-head saturation-test comparison on paired simulations ({model_label(model)})</text>',
    ]
    panel_w, panel_h = 330, 235
    x0, y0 = 76, 94
    x_gap, y_gap = 42, 112
    for col, (formula, formula_label) in enumerate(FORMULAS):
        svg.append(f'<text x="{x0 + col * (panel_w + x_gap) + panel_w / 2:.1f}" y="61" text-anchor="middle" font-size="14" font-weight="600">{formula_label}</text>')
        for row_idx, (tree_case, tree_label) in enumerate(tree_cases):
            panel(
                svg,
                rows,
                formula,
                tree_case,
                tree_label,
                x0 + col * (panel_w + x_gap),
                y0 + row_idx * (panel_h + y_gap),
                panel_w,
                panel_h,
                show_y_label=(col == 0),
                show_x_label=(row_idx == len(tree_cases) - 1),
            )
    legend(svg, 1200, 110)
    note = "Full paper grid, paired by alignment, target branch, scenario and threshold."
    if problems:
        note = "Incomplete plotting-development output; not for manuscript use."
    svg.append(f'<text x="76" y="{height - 30}" font-size="11" fill="#555">{note}</text>')
    svg.append("</svg>")
    output.write_text("\n".join(svg) + "\n", encoding="utf-8")


def maybe_convert_pdf(svg_path):
    pdf_path = svg_path.with_suffix(".pdf")
    try:
        subprocess.run(["rsvg-convert", "-f", "pdf", "-o", str(pdf_path), str(svg_path)], check=True)
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None
    return pdf_path


def main():
    parser = argparse.ArgumentParser(description="Plot full head-to-head SatuTe simulation curves.")
    parser.add_argument("--summary", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--expected-reps", type=int, default=1000)
    parser.add_argument("--tree-cases", default="five_external,sixteen_internal")
    parser.add_argument("--model", default="", help="Model to plot. Use JC or GTR_PF06346 for the bundled examples.")
    parser.add_argument("--allow-incomplete", action="store_true", help="Only for development plots; manuscript plots require complete data.")
    args = parser.parse_args()

    rows = load_rows(args.summary)
    model = resolve_model(args.model, rows)
    rows = [row for row in rows if row["model"] == model]
    tree_cases = selected_tree_cases(args.tree_cases)
    problems = validate_full(rows, args.expected_reps, args.allow_incomplete, tree_cases)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    suffix = model_label(model)
    baseline = outdir / f"figure_simulated_data_main_reproduction_{suffix}.svg"
    comparison = outdir / f"figure_head_to_head_formula_comparison_{suffix}.svg"
    write_baseline_figure(rows, baseline, problems, tree_cases, model)
    write_formula_grid(rows, comparison, problems, tree_cases, model)
    baseline_pdf = maybe_convert_pdf(baseline)
    comparison_pdf = maybe_convert_pdf(comparison)

    print(f"Baseline SVG:   {baseline}")
    if baseline_pdf:
        print(f"Baseline PDF:   {baseline_pdf}")
    print(f"Comparison SVG: {comparison}")
    if comparison_pdf:
        print(f"Comparison PDF: {comparison_pdf}")
    if problems:
        print(f"WARNING: incomplete development plot with {len(problems)} validation problem(s).")


if __name__ == "__main__":
    main()
