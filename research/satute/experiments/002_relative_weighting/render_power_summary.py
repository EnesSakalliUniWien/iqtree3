#!/usr/bin/env python3

import argparse
import csv
import math
import sys
from pathlib import Path
from xml.sax.saxutils import escape


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from satute_analysis.publication_figures import (  # noqa: E402
    CROSS_MODEL_CONTRACT,
    stable_case_summary,
    stable_figure_root,
)
from satute_analysis.benchmark_contracts import applicable_decision_rules  # noqa: E402
from tools.figures.export import export_svg_pdf  # noqa: E402

OUTDIR = stable_figure_root(CROSS_MODEL_CONTRACT.id)

DATASETS = [
    {
        "label": "JC",
        "path": stable_case_summary("jc_JC"),
        "simulation_model": "JC",
        "evaluation_model": "JC",
        "site_lengths": [str(value) for value in CROSS_MODEL_CONTRACT.site_lengths],
        "branch_ticks": [0.1, 0.5, 1.0, 2.5, 5.0, 10.0],
    },
    {
        "label": "GTR+F+G4",
        "path": stable_case_summary("gtr_GTR_PF06346_G4"),
        "simulation_model": "GTR{0.6676,3.7807,4.2833,0.5354,0.8718,1.0}+F{0.125,0.436,0.191,0.245}+G4{0.5}",
        "evaluation_model": "GTR{0.6676,3.7807,4.2833,0.5354,0.8718,1.0}+F{0.125,0.436,0.191,0.245}+G4{0.5}",
        "site_lengths": [str(value) for value in CROSS_MODEL_CONTRACT.site_lengths],
        "branch_ticks": [4, 6, 8, 10, 12],
    },
    {
        "label": "LG+G4",
        "path": stable_case_summary("protein_LG_G4"),
        "simulation_model": "LG+G4{0.5}",
        "evaluation_model": "LG+G4{0.5}",
        "site_lengths": [str(value) for value in CROSS_MODEL_CONTRACT.site_lengths],
        "branch_ticks": [4, 6, 8, 10, 12],
    },
]

TREE_CASES = [
    ("five_external", "five-taxon external branch"),
    ("sixteen_internal", "16-taxon internal branch"),
]

FORMULAS = [
    ("dominant", "Published SatuTe", "#0A3FBA", "circle"),
    ("eigenvalue_weighted", "Eigenvalue-weighted SatuTe", "#E69504", "diamond"),
]

DECISION_RULE_LABELS = {
    "unadjusted": "unadjusted",
    "taxon_bonferroni": "taxon-pair Bonferroni",
    "by_fdr": "Benjamini--Yekutieli FDR",
}

SITE_DASH = {
    "100": "8,5",
    "250": "4,4",
    "1000": "",
    "10000": "1,5",
}


def read_rows(path, decision_rule):
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SystemExit(f"Summary is empty: {path}")
    if "decision_rule" not in rows[0]:
        raise SystemExit(
            f"Summary must use schema version 2 with a decision_rule column: {path}"
        )
    rows = [row for row in rows if row["decision_rule"] == decision_rule]
    if not rows:
        raise SystemExit(f"Summary contains no rows for decision rule {decision_rule}: {path}")
    return rows


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


def grouped(rows, dataset, tree_case, scenario):
    values = {}
    for row in rows:
        if row["simulation_model"] != dataset["simulation_model"]:
            continue
        if row["evaluation_model"] != dataset["evaluation_model"]:
            continue
        if row["scenario"] != scenario:
            continue
        if row["tree_case"] != tree_case:
            continue
        if row["fraction_informative"] == "":
            continue
        key = (row["formula"], row["nsites"])
        branch = float(row["branch_length"])
        series = values.setdefault(key, {})
        if branch in series:
            raise SystemExit(
                "Duplicate summary cell after decision-rule filtering: "
                f"{dataset['label']}/{tree_case}/{row['formula']}/n{row['nsites']}/b{branch}"
            )
        series[branch] = float(row["fraction_informative"])
    return values


def draw_panel(
    svg,
    rows,
    dataset,
    tree_case,
    title,
    left,
    top,
    width,
    height,
    show_y,
    show_x,
    scenario,
):
    data = grouped(rows, dataset, tree_case, scenario)
    all_branch_lengths = sorted({branch for series in data.values() for branch in series})
    if not all_branch_lengths:
        raise SystemExit(
            f"No {scenario} rows for {dataset['label']}/{tree_case}"
        )
    xmin, xmax = min(all_branch_lengths), max(all_branch_lengths)
    if xmin == xmax:
        center = math.log10(xmin)
        xmin, xmax = 10 ** (center - 0.5), 10 ** (center + 0.5)

    svg.append(f'<text x="{left + width / 2:.1f}" y="{top - 18}" text-anchor="middle" font-size="17" font-weight="700">{escape(title)}</text>')
    svg.append(f'<rect x="{left}" y="{top}" width="{width}" height="{height}" fill="white" stroke="#444" stroke-width="1.1"/>')

    for value in [0, 0.25, 0.5, 0.75, 1.0]:
        y = sy(value, top, height)
        svg.append(f'<line x1="{left}" x2="{left + width}" y1="{y:.2f}" y2="{y:.2f}" stroke="#e6e6e6" stroke-width="1"/>')
        if show_y:
            svg.append(f'<text x="{left - 10}" y="{y + 5:.2f}" text-anchor="end" font-size="14">{value:g}</text>')

    for tick in dataset["branch_ticks"]:
        if tick < xmin or tick > xmax:
            continue
        x = sx(tick, xmin, xmax, left, width)
        svg.append(f'<line x1="{x:.2f}" x2="{x:.2f}" y1="{top}" y2="{top + height}" stroke="#eeeeee" stroke-width="1"/>')
        if show_x:
            svg.append(f'<text x="{x:.2f}" y="{top + height + 22}" text-anchor="middle" font-size="14">{tick:g}</text>')

    for formula, _label, color, marker in FORMULAS:
        for nsites in dataset["site_lengths"]:
            series = data.get((formula, nsites), {})
            points = [(sx(branch, xmin, xmax, left, width), sy(series[branch], top, height)) for branch in all_branch_lengths if branch in series]
            if not points:
                continue
            dash = SITE_DASH.get(nsites, "")
            dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
            stroke_width = 4.6 if formula == "dominant" else 2.4
            svg.append(
                f'<path d="{path_data(points)}" fill="none" stroke="{color}" '
                f'stroke-width="{stroke_width}"{dash_attr}/>'
            )
            for x, y in points:
                if marker == "circle":
                    svg.append(
                        f'<circle cx="{x:.2f}" cy="{y:.2f}" r="4.1" fill="white" '
                        f'stroke="{color}" stroke-width="2.0"/>'
                    )
                else:
                    size = 3.7
                    svg.append(
                        f'<polygon points="{x:.2f},{y-size:.2f} {x+size:.2f},{y:.2f} '
                        f'{x:.2f},{y+size:.2f} {x-size:.2f},{y:.2f}" fill="white" '
                        f'stroke="{color}" stroke-width="1.8"/>'
                    )

    deltas = []
    for nsites in dataset["site_lengths"]:
        published = data.get(("dominant", nsites), {})
        weighted = data.get(("eigenvalue_weighted", nsites), {})
        for branch in set(published) & set(weighted):
            deltas.append(weighted[branch] - published[branch])
    maximum_delta = max(deltas) if deltas else float("nan")
    annotation = (
        "curves coincide (Δ = 0)"
        if math.isclose(maximum_delta, 0.0, abs_tol=1e-12)
        else f"max Δ = {maximum_delta:.3f}"
    )
    svg.append(
        f'<text x="{left + width - 10}" y="{top + 20}" text-anchor="end" '
        f'font-size="13" font-weight="700" fill="#444444">{annotation}</text>'
    )

    if show_y:
        x = left - 58
        y = top + height / 2
        svg.append(f'<text x="{x}" y="{y:.1f}" text-anchor="middle" font-size="15" transform="rotate(-90 {x} {y:.1f})">fraction phylogenetically informative</text>')
    if show_x:
        svg.append(f'<text x="{left + width / 2:.1f}" y="{top + height + 50}" text-anchor="middle" font-size="15">target branch length</text>')


def legend(svg, x, y, site_lengths):
    svg.append(f'<text x="{x}" y="{y}" font-size="17" font-weight="700">Formula</text>')
    for i, (_formula, label, color, marker) in enumerate(FORMULAS):
        yy = y + 28 + i * 28
        width = 4.6 if marker == "circle" else 2.4
        svg.append(
            f'<line x1="{x}" x2="{x + 42}" y1="{yy}" y2="{yy}" '
            f'stroke="{color}" stroke-width="{width}"/>'
        )
        if marker == "circle":
            svg.append(
                f'<circle cx="{x + 21}" cy="{yy}" r="4.1" fill="white" '
                f'stroke="{color}" stroke-width="2"/>'
            )
        else:
            svg.append(
                f'<polygon points="{x+21},{yy-4} {x+25},{yy} {x+21},{yy+4} '
                f'{x+17},{yy}" fill="white" stroke="{color}" stroke-width="1.8"/>'
            )
        svg.append(f'<text x="{x + 54}" y="{yy + 5}" font-size="15">{escape(label)}</text>')

    y2 = y + 100
    svg.append(f'<text x="{x}" y="{y2}" font-size="17" font-weight="700">Alignment length</text>')
    for i, nsites in enumerate(site_lengths):
        yy = y2 + 28 + i * 28
        dash = SITE_DASH.get(nsites, "")
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        svg.append(f'<line x1="{x}" x2="{x + 42}" y1="{yy}" y2="{yy}" stroke="#222222" stroke-width="2.8"{dash_attr}/>')
        svg.append(f'<text x="{x + 54}" y="{yy + 5}" font-size="15">{nsites} sites</text>')


def write_recalculation_tables(
    datasets, rows_by_path, scenario, decision_rule, output_dir
):
    selected_rows = []
    claim_rows = []
    for dataset in datasets:
        rows = [
            row
            for row in rows_by_path[dataset["path"]]
            if row["simulation_model"] == dataset["simulation_model"]
            and row["evaluation_model"] == dataset["evaluation_model"]
            and row["scenario"] == scenario
            and row["decision_rule"] == decision_rule
            and row["nsites"] in dataset["site_lengths"]
            and row["tree_case"] in {value for value, _label in TREE_CASES}
        ]
        for row in rows:
            selected_rows.append({"model_label": dataset["label"], **row})
        by_cell = {}
        for row in rows:
            key = (row["tree_case"], row["nsites"], row["branch_length"])
            by_cell.setdefault(key, {})[row["formula"]] = float(
                row["fraction_informative"]
            )
        for tree_case, _tree_label in TREE_CASES:
            candidates = []
            for (cell_tree, nsites, branch_length), formulas in by_cell.items():
                if cell_tree != tree_case or set(formulas) != {
                    "dominant",
                    "eigenvalue_weighted",
                }:
                    continue
                candidates.append(
                    (
                        formulas["eigenvalue_weighted"] - formulas["dominant"],
                        int(nsites),
                        float(branch_length),
                        formulas["dominant"],
                        formulas["eigenvalue_weighted"],
                    )
                )
            if not candidates:
                raise SystemExit(
                    f"No paired Figure 6 cells for {dataset['label']}/{tree_case}"
                )
            delta, nsites, branch_length, published, weighted = max(candidates)
            claim_rows.append(
                {
                    "model": dataset["label"],
                    "tree_case": tree_case,
                    "scenario": scenario,
                    "decision_rule": decision_rule,
                    "nsites": nsites,
                    "branch_length": branch_length,
                    "published_fraction": published,
                    "eigenvalue_weighted_fraction": weighted,
                    "delta_fraction": delta,
                }
            )

    data_fields = [
        "model_label",
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
    ]
    with (output_dir / "figure_power_curve_summary_data.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=data_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(
            {field: row[field] for field in data_fields} for row in selected_rows
        )
    claim_fields = list(claim_rows[0])
    with (output_dir / "figure_power_curve_summary_claims.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=claim_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(claim_rows)


def main():
    parser = argparse.ArgumentParser(description="Render the reviewed relative-weighting power summary.")
    parser.add_argument("--output-dir", type=Path, default=OUTDIR)
    parser.add_argument("--jc-summary", type=Path, default=DATASETS[0]["path"])
    parser.add_argument(
        "--gtr-summary",
        "--rate-summary",
        dest="gtr_summary",
        type=Path,
        default=DATASETS[1]["path"],
        help="Summary for the GTR+F+G4 dataset (legacy alias: --rate-summary).",
    )
    parser.add_argument("--lg-summary", type=Path, default=DATASETS[2]["path"])
    parser.add_argument("--gtr-simulation-model", default=DATASETS[1]["simulation_model"])
    parser.add_argument("--gtr-evaluation-model", default=DATASETS[1]["evaluation_model"])
    parser.add_argument("--lg-simulation-model", default=DATASETS[2]["simulation_model"])
    parser.add_argument("--lg-evaluation-model", default=DATASETS[2]["evaluation_model"])
    parser.add_argument(
        "--decision-rule",
        choices=tuple(DECISION_RULE_LABELS),
        default=CROSS_MODEL_CONTRACT.decision_rule,
        help="Decision contract to render; schema-v2 summaries contain all rules.",
    )
    parser.add_argument(
        "--scenario",
        choices=("true_tree_fixed_lengths", "true_topology_ml_lengths", "ml_tree"),
        default=CROSS_MODEL_CONTRACT.scenario,
    )
    parser.add_argument("--expected-reps", type=int, default=1000)
    parser.add_argument(
        "--allow-incomplete",
        action="store_true",
        help="Allow development-only summaries with incomplete replicate cells.",
    )
    args = parser.parse_args()
    allowed_rules = {rule.name for rule in applicable_decision_rules(args.scenario)}
    if args.decision_rule not in allowed_rules:
        raise SystemExit(
            f"{args.decision_rule} is not applicable to scenario {args.scenario}"
        )
    datasets = [dict(dataset) for dataset in DATASETS]
    datasets[0]["path"] = args.jc_summary
    datasets[1]["path"] = args.gtr_summary
    datasets[2]["path"] = args.lg_summary
    datasets[1]["simulation_model"] = args.gtr_simulation_model
    datasets[1]["evaluation_model"] = args.gtr_evaluation_model
    datasets[2]["simulation_model"] = args.lg_simulation_model
    datasets[2]["evaluation_model"] = args.lg_evaluation_model
    figure_dir = args.output_dir / "figures"
    data_dir = args.output_dir / "data"
    out_svg = figure_dir / "figure_power_curve_summary.svg"
    out_pdf = figure_dir / "figure_power_curve_summary.pdf"

    rows_by_path = {}
    for dataset in datasets:
        rows_by_path.setdefault(
            dataset["path"], read_rows(dataset["path"], args.decision_rule)
        )

    if args.expected_reps < 1:
        raise SystemExit("--expected-reps must be positive")
    for path, rows in rows_by_path.items():
        problems = []
        for row in rows:
            try:
                total = int(row["evaluated"]) + int(row["missing_split"])
            except (KeyError, TypeError, ValueError) as exc:
                raise SystemExit(f"Malformed summary row in {path}: {row}") from exc
            if total != args.expected_reps:
                problems.append(
                    f"{row.get('tree_case')}/n{row.get('nsites')}/b{row.get('branch_length')} "
                    f"has evaluated+missing_split={total}, expected {args.expected_reps}"
                )
        if problems and not args.allow_incomplete:
            sample = "\n".join(f"- {problem}" for problem in problems[:10])
            raise SystemExit(
                f"Summary {path} is not complete for {args.expected_reps} replicates. "
                "Use --allow-incomplete only for development plots.\n" + sample
            )

    site_lengths = {dataset["label"]: dataset["site_lengths"] for dataset in datasets}
    tree_case_sets = [
        {row["tree_case"] for row in rows_by_path[dataset["path"]]}
        for dataset in datasets
    ]
    tree_cases = [
        ("five_external", "five-taxon external branch"),
        ("sixteen_internal", "16-taxon internal branch"),
    ]
    tree_cases = [
        item for item in tree_cases if all(item[0] in available for available in tree_case_sets)
    ]
    if not tree_cases:
        raise SystemExit("No tree case is present in every selected summary")

    figure_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    width, height = 1480, 1300
    panel_w, panel_h = 395, 220
    x0, y0 = 135, 130
    xgap, ygap = 105, 170

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{width / 2:.1f}" y="46" text-anchor="middle" font-size="30" font-weight="700">Cross-model SatuTe power under inferred-tree testing</text>',
        f'<text x="{width / 2:.1f}" y="78" text-anchor="middle" font-size="17" fill="#333333">ML-tree scenario; {args.expected_reps:,} paired replicates per point; {DECISION_RULE_LABELS[args.decision_rule]}</text>',
    ]

    for row_idx, dataset in enumerate(datasets):
        top = y0 + row_idx * (panel_h + ygap)
        svg.append(f'<text x="30" y="{top + panel_h / 2:.1f}" text-anchor="middle" font-size="21" font-weight="700" transform="rotate(-90 30 {top + panel_h / 2:.1f})">{escape(dataset["label"])}</text>')
        rows = rows_by_path[dataset["path"]]
        for col_idx, (tree_case, tree_label) in enumerate(tree_cases):
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
                scenario=args.scenario,
            )

    common_site_lengths = sorted(
        {site for values in site_lengths.values() for site in values}, key=int
    )
    legend(svg, 1160, 140, common_site_lengths)
    svg.append(
        f'<text x="135" y="1250" font-size="14" fill="#555555">'
        f'Only the common 100- and 1,000-site grids are shown. '
        f'Curves use {escape(DECISION_RULE_LABELS[args.decision_rule])}.</text>'
    )
    svg.append("</svg>")

    out_svg.write_text("\n".join(svg) + "\n", encoding="utf-8")
    export_svg_pdf(out_svg, out_pdf)
    write_recalculation_tables(
        datasets,
        rows_by_path,
        args.scenario,
        args.decision_rule,
        data_dir,
    )
    print(out_svg)
    print(out_pdf)


if __name__ == "__main__":
    main()
