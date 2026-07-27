#!/usr/bin/env python3
"""Validate, summarize and render the HCV codon-position window example."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path

from render_curated import convert_svg, line, polyline, rect, svg_document, text


BRANCHES = {
    tuple(sorted(("HCV_1a_H77", "HCV_1b_Con1"))): "genotype 1",
    tuple(sorted(("HCV_2a_JFH1", "HCV_2b_HCJ8"))): "genotype 2",
    tuple(sorted(("HCV_3a_NZL1", "HCV_3b_TrKj"))): "genotype 3",
}
BRANCH_ORDER = ("genotype 1", "genotype 2", "genotype 3")
FORMULAS = ("dominant", "eigenvalue_weighted")
FORMULA_LABELS = {
    "dominant": "published SatuTe",
    "eigenvalue_weighted": "eigenvalue-weighted SatuTe",
}
COLORS = {"dominant": "#3569B0", "eigenvalue_weighted": "#D67B16"}
PARTITION_LABELS = {"pos12": "codon positions 1+2", "pos3": "codon position 3"}
GENES = (
    ("Core", 1, 191),
    ("E1", 192, 383),
    ("E2", 384, 746),
    ("p7", 747, 809),
    ("NS2", 810, 1026),
    ("NS3", 1027, 1657),
    ("NS4A", 1658, 1711),
    ("NS4B", 1712, 1972),
    ("NS5A", 1973, 2420),
    ("NS5B", 2421, 3011),
)
Z_ALPHA_001 = 2.3263478740408408


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pos12-detail", type=Path)
    parser.add_argument("--pos3-detail", type=Path)
    parser.add_argument("--window-manifest", type=Path)
    parser.add_argument(
        "--validated-detail",
        type=Path,
        help="Re-render a previously validated hcv_window_detail.tsv release table.",
    )
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    raw_inputs = (args.pos12_detail, args.pos3_detail, args.window_manifest)
    if args.validated_detail:
        if any(raw_inputs):
            parser.error("--validated-detail cannot be combined with raw detail inputs")
    elif not all(raw_inputs):
        parser.error(
            "provide --validated-detail, or provide --pos12-detail, "
            "--pos3-detail and --window-manifest"
        )
    return args


def gene_for(position: int) -> str:
    for gene, start, end in GENES:
        if start <= position <= end:
            return gene
    return "alignment gap"


def read_manifest(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    manifest = {Path(row["filename"]).stem: row for row in rows}
    if len(manifest) != len(rows):
        raise ValueError("Window manifest contains duplicate filenames")
    counts = Counter(row["partition"] for row in rows)
    if counts != Counter({"pos12": 507, "pos3": 252}):
        raise ValueError(f"Unexpected HCV window grid: {counts}")
    return manifest


def selected_rows(detail_paths: dict[str, Path], manifest: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    output: list[dict[str, str]] = []
    observed = Counter()
    for partition, path in detail_paths.items():
        with path.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                split = tuple(sorted(part.strip() for part in row["split"].split(",") if part.strip()))
                if split not in BRANCHES:
                    continue
                if row["formula"] not in FORMULAS:
                    continue
                window = manifest.get(row["dataset"])
                if window is None:
                    raise ValueError(f"Dataset absent from window manifest: {row['dataset']}")
                if window["partition"] != partition:
                    raise ValueError(f"Partition mismatch for {row['dataset']}")
                midpoint = int(window["h77_polyprotein_midpoint"])
                selected = {
                    "partition": partition,
                    "dataset": row["dataset"],
                    "start_site": window["start_site"],
                    "end_site": window["end_site"],
                    "h77_polyprotein_midpoint": str(midpoint),
                    "gene": gene_for(midpoint),
                    "branch": BRANCHES[split],
                    "formula": row["formula"],
                    "branch_length": row["branch_length"],
                    "satZ": row["satZ"],
                    "satP": row["satP"],
                    "decision": row["decision"],
                    "decision_taxon_bonf": row["decision_taxon_bonf"],
                }
                output.append(selected)
                observed[(partition, BRANCHES[split], row["formula"])] += 1

    expected = Counter()
    for partition, count in (("pos12", 507), ("pos3", 252)):
        for branch in BRANCH_ORDER:
            for formula in FORMULAS:
                expected[(partition, branch, formula)] = count
    if observed != expected:
        raise ValueError(f"Incomplete HCV branch/formula/window grid: {observed - expected}; {expected - observed}")
    return sorted(
        output,
        key=lambda row: (
            row["partition"],
            BRANCH_ORDER.index(row["branch"]),
            FORMULAS.index(row["formula"]),
            int(row["start_site"]),
        ),
    )


def read_validated_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "partition",
        "dataset",
        "start_site",
        "end_site",
        "h77_polyprotein_midpoint",
        "gene",
        "branch",
        "formula",
        "branch_length",
        "satZ",
        "satP",
        "decision",
        "decision_taxon_bonf",
    }
    missing = required - set(rows[0] if rows else ())
    if missing:
        raise ValueError(f"Validated HCV detail is missing columns: {sorted(missing)}")
    observed = Counter(
        (row["partition"], row["branch"], row["formula"]) for row in rows
    )
    expected = Counter()
    for partition, count in (("pos12", 507), ("pos3", 252)):
        for branch in BRANCH_ORDER:
            for formula in FORMULAS:
                expected[(partition, branch, formula)] = count
    if observed != expected:
        raise ValueError(
            f"Incomplete validated HCV grid: {observed - expected}; {expected - observed}"
        )
    return sorted(
        rows,
        key=lambda row: (
            row["partition"],
            BRANCH_ORDER.index(row["branch"]),
            FORMULAS.index(row["formula"]),
            int(row["start_site"]),
        ),
    )


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows: list[dict[str, str]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    genes: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["partition"], row["branch"], row["formula"])].append(row)
        genes[(row["partition"], row["gene"], row["branch"], row["formula"])].append(row)

    summary: list[dict[str, object]] = []
    for (partition, branch, formula), selected in sorted(grouped.items()):
        zscores = [float(row["satZ"]) for row in selected]
        saturated = sum(float(row["satZ"]) < Z_ALPHA_001 for row in selected)
        summary.append(
            {
                "partition": partition,
                "branch": branch,
                "formula": formula,
                "windows": len(selected),
                "saturated_windows": saturated,
                "saturated_fraction": saturated / len(selected),
                "min_satZ": min(zscores),
                "median_satZ": statistics.median(zscores),
                "max_satZ": max(zscores),
                "alpha": 0.01,
                "z_threshold": Z_ALPHA_001,
            }
        )

    gene_summary: list[dict[str, object]] = []
    for (partition, gene, branch, formula), selected in sorted(genes.items()):
        saturated = sum(float(row["satZ"]) < Z_ALPHA_001 for row in selected)
        gene_summary.append(
            {
                "partition": partition,
                "gene": gene,
                "branch": branch,
                "formula": formula,
                "windows": len(selected),
                "saturated_windows": saturated,
                "saturated_fraction": saturated / len(selected),
            }
        )
    return summary, gene_summary


def render_curves(rows: list[dict[str, str]], output_base: Path) -> None:
    by_key: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_key[(row["partition"], row["branch"], row["formula"])].append(row)

    width, height = 1800, 1420
    left_margin, right_margin = 105, 35
    column_gap = 90
    panel_width = (width - left_margin - right_margin - column_gap) / 2
    panel_height = 300
    row_tops = (130, 555, 980)
    y_min, y_max = -0.5, 10.8
    body: list[str] = []

    for index, formula in enumerate(FORMULAS):
        x = 630 + index * 330
        body.append(line(x, 48, x + 70, 48, color=COLORS[formula], width=4.0))
        body.append(text(x + 84, 55, FORMULA_LABELS[formula], size=22, anchor="start"))
    body.append(line(1280, 48, 1350, 48, color="#666666", width=2.0, dash="10 7"))
    body.append(text(1364, 55, "alpha = 0.01", size=22, anchor="start"))

    for column, partition in enumerate(("pos12", "pos3")):
        left = left_margin + column * (panel_width + column_gap)
        body.append(text(left + panel_width / 2, 94, PARTITION_LABELS[partition], size=27, weight="700"))
        for row_index, branch in enumerate(BRANCH_ORDER):
            top = row_tops[row_index]
            bottom = top + panel_height
            body.append(rect(left, top, panel_width, panel_height, fill="#ffffff", stroke="#222222"))

            def sx(value: float) -> float:
                return left + panel_width * (value - 1) / 3010

            def sy(value: float) -> float:
                clipped = min(y_max, max(y_min, value))
                return bottom - panel_height * (clipped - y_min) / (y_max - y_min)

            for gene_index, (gene, start, end) in enumerate(GENES):
                x1, x2 = sx(start), sx(end)
                if gene_index % 2:
                    body.append(rect(x1, top, x2 - x1, panel_height, fill="#F4F6F8"))
                body.append(line(x1, top, x1, bottom, color="#D5D9DE", width=0.8))
                if row_index == len(BRANCH_ORDER) - 1 and end - start >= 110:
                    body.append(text((x1 + x2) / 2, bottom + 28, gene, size=15))

            for tick in (0, 2, 4, 6, 8, 10):
                y = sy(tick)
                body.append(line(left, y, left + panel_width, y, color="#D8DCE0", width=0.8))
                if column == 0:
                    body.append(text(left - 13, y + 5, str(tick), size=15, anchor="end"))
            body.append(line(left, sy(Z_ALPHA_001), left + panel_width, sy(Z_ALPHA_001), color="#666666", width=1.7, dash="10 7"))
            for formula in FORMULAS:
                selected = sorted(
                    by_key[(partition, branch, formula)], key=lambda row: int(row["start_site"])
                )
                points = [
                    (sx(float(row["h77_polyprotein_midpoint"])), sy(float(row["satZ"])))
                    for row in selected
                ]
                body.append(polyline(points, color=COLORS[formula], width=2.0))

            if column == 0:
                body.append(text(34, top + panel_height / 2, branch, size=21, weight="700", rotate=-90))
            if row_index == 1 and column == 0:
                body.append(text(78, top + panel_height / 2, "Window z-score", size=18, rotate=-90))

    body.append(text(width / 2, height - 20, "H77 polyprotein residue", size=22))
    svg = svg_document(width, height, body)
    output_base.with_suffix(".svg").write_text(svg, encoding="utf-8")
    convert_svg(svg, output_base, 2700)


def render_bars(summary: list[dict[str, object]], output_base: Path) -> None:
    by_key = {
        (str(row["partition"]), str(row["branch"]), str(row["formula"])): row
        for row in summary
    }
    width, height = 1450, 760
    left, right, top, bottom = 115, 35, 105, 150
    plot_width, plot_height = width - left - right, height - top - bottom
    y_max = 16.0
    body: list[str] = [rect(left, top, plot_width, plot_height, fill="#ffffff", stroke="#222222")]
    for tick in (0, 4, 8, 12, 16):
        y = top + plot_height * (1 - tick / y_max)
        body.append(line(left, y, left + plot_width, y, color="#D8DCE0", width=0.9))
        body.append(text(left - 14, y + 5, str(tick), size=16, anchor="end"))

    group_width = plot_width / len(BRANCH_ORDER)
    bar_width = 72
    offsets = (-126, -42, 42, 126)
    combinations = (("pos12", "dominant"), ("pos12", "eigenvalue_weighted"), ("pos3", "dominant"), ("pos3", "eigenvalue_weighted"))
    fills = ("#ADC5E6", "#F1C28C", COLORS["dominant"], COLORS["eigenvalue_weighted"])
    for branch_index, branch in enumerate(BRANCH_ORDER):
        center = left + group_width * (branch_index + 0.5)
        for offset, (partition, formula), fill in zip(offsets, combinations, fills, strict=True):
            row = by_key[(partition, branch, formula)]
            value = 100 * float(row["saturated_fraction"])
            x = center + offset - bar_width / 2
            height_value = plot_height * value / y_max
            body.append(rect(x, top + plot_height - height_value, bar_width, height_value, fill=fill))
            body.append(text(x + bar_width / 2, top + plot_height - height_value - 9, f"{value:.1f}", size=15))
        body.append(text(center, top + plot_height + 35, branch, size=20, weight="700"))

    legend = (
        ("positions 1+2, published SatuTe", "#ADC5E6"),
        ("positions 1+2, eigenvalue-weighted", "#F1C28C"),
        ("position 3, published SatuTe", COLORS["dominant"]),
        ("position 3, eigenvalue-weighted", COLORS["eigenvalue_weighted"]),
    )
    for index, (label, fill) in enumerate(legend):
        x = 210 + index * 305
        body.append(rect(x, 43, 28, 22, fill=fill))
        body.append(text(x + 40, 61, label, size=16, anchor="start"))
    body.append(text(35, top + plot_height / 2, "Saturated windows (%)", size=20, rotate=-90))
    body.append(text(width / 2, height - 22, "Equal 36-nucleotide windows; alpha = 0.01", size=18))
    svg = svg_document(width, height, body)
    output_base.with_suffix(".svg").write_text(svg, encoding="utf-8")
    convert_svg(svg, output_base, 2175)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    if args.validated_detail:
        rows = read_validated_rows(args.validated_detail)
    else:
        manifest = read_manifest(args.window_manifest)
        rows = selected_rows(
            {"pos12": args.pos12_detail, "pos3": args.pos3_detail}, manifest
        )
    summary, gene_summary = summarize(rows)
    write_tsv(args.outdir / "hcv_window_detail.tsv", rows, list(rows[0]))
    write_tsv(args.outdir / "hcv_window_summary.tsv", summary, list(summary[0]))
    write_tsv(args.outdir / "hcv_window_gene_summary.tsv", gene_summary, list(gene_summary[0]))
    render_curves(rows, args.outdir / "hcv_codon_position_window_curves")
    render_bars(summary, args.outdir / "hcv_codon_position_saturated_windows")
    print(f"Wrote and validated {len(rows)} branch/formula/window rows in {args.outdir}")


if __name__ == "__main__":
    main()
