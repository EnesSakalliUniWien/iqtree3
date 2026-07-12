#!/usr/bin/env python3
"""Run and plot focused diagnostics for the three native SatuTe formulas.

The experiment fixes the model and tree so that movement in satC, satSE, satZ,
and the decision rate can be attributed to the statistic rather than parameter
estimation.  A saturated-branch sample is split into threshold-training and
held-out calibration halves.
"""

import argparse
import csv
import hashlib
import math
import os
import shutil
import statistics
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from xml.sax.saxutils import escape


PROJECT_ROOT = Path(__file__).resolve().parents[2]
IQTREE_DEFAULT = PROJECT_ROOT.parents[1] / "build" / "iqtree3"


DEFAULT_MODEL = "GTR{1,2,1,1,2,1}+F{0.30,0.20,0.20,0.30}+G4{0.5}"
DEFAULT_BRANCHES = "0.1,0.5,1,2,5,10,20,30,50,75,100,150,200,500,10000"
FORMULAS = (
    "dominant",
    "eigenvalue_weighted",
    "mixture_likelihood_weighted",
)
FORMULA_STYLE = {
    "dominant": ("Dominant", "#0A3FBA", "", "circle"),
    "eigenvalue_weighted": ("Relative eigenvalue-weighted", "#E69504", "9,5", "square"),
    "mixture_likelihood_weighted": ("Soft likelihood mixture", "#B6497D", "2,4", "triangle"),
}
DETAIL_FIELDS = (
    "branch_length",
    "replicate",
    "seed",
    "formula",
    "satC",
    "satVar",
    "satSE",
    "satZ",
    "satP",
    "information_fraction",
    "saturation_index",
    "decision",
    "valid_sites",
    "skipped_sites",
    "modes",
    "eigenvalues",
    "weights",
)
CATEGORY_FIELDS = (
    "branch_length",
    "replicate",
    "seed",
    "formula",
    "category",
    "rate_multiplier",
    "rate_sites",
    "valid_sites",
    "skipped_sites",
    "satC",
    "satVar",
    "satSE",
    "satZ",
    "weights",
)


def parse_numbers(text, cast=float):
    return [cast(value.strip()) for value in text.split(",") if value.strip()]


def branch_tag(value):
    return f"{value:.10g}".replace("-", "m").replace(".", "p")


def simulation_seed(base_seed, sites, branch_index, replicate):
    """Return a reproducible, nonconsecutive positive 31-bit simulator seed."""
    payload = f"{base_seed}:{sites}:{branch_index}:{replicate}".encode("ascii")
    value = int.from_bytes(hashlib.blake2s(payload, digest_size=8).digest(), "big")
    return 1 + value % 2147483646


def run_command(command):
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode:
        rendered = " ".join(command)
        raise RuntimeError(f"Command failed ({result.returncode}): {rendered}\n{result.stderr[-4000:]}")


def parse_stat(path, target_taxa=("A", "B")):
    target = set(target_taxa)
    selected = []
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader((line for line in handle if line.strip() and not line.startswith("#")), delimiter="\t")
        for row in reader:
            split = {item for item in row.get("Split", "").split(",") if item}
            if split == target:
                selected.append(row)
    pooled = {row["Formula"]: row for row in selected if row["RateCategory"] == "pooled"}
    missing = set(FORMULAS) - set(pooled)
    if missing:
        raise ValueError(f"Missing pooled formulas in {path}: {sorted(missing)}")
    return pooled, selected


def parse_rate_categories(path):
    categories = []
    in_table = False
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            fields = raw.split()
            if fields[:3] == ["Category", "Relative_rate", "Proportion"]:
                in_table = True
                continue
            if not in_table:
                continue
            if len(fields) == 3 and fields[0].isdigit():
                try:
                    categories.append(
                        {
                            "category": fields[0],
                            "rate": float(fields[1]),
                            "proportion": float(fields[2]),
                        }
                    )
                    continue
                except ValueError:
                    pass
            if categories:
                break
    if not categories:
        raise ValueError(f"Could not parse the rate-category table from {path}")
    total = sum(item["proportion"] for item in categories)
    for item in categories:
        item["proportion"] /= total
    return categories


def pooled_to_detail(branch_length, replicate, seed, formula, row):
    return {
        "branch_length": f"{branch_length:.10g}",
        "replicate": str(replicate),
        "seed": str(seed),
        "formula": formula,
        "satC": row["satC"],
        "satVar": row["satVar"],
        "satSE": row["satSE"],
        "satZ": row["satZ"],
        "satP": row["satP"],
        "information_fraction": row["InformationFraction"],
        "saturation_index": row["SaturationIndex"],
        "decision": row["Decision"],
        "valid_sites": row["ValidSites"],
        "skipped_sites": row["SkippedSites"],
        "modes": row["Modes"],
        "eigenvalues": row["Eigenvalues"],
        "weights": row["Weights"],
    }


def native_category_to_detail(branch_length, replicate, seed, row):
    return {
        "branch_length": f"{branch_length:.10g}",
        "replicate": str(replicate),
        "seed": str(seed),
        "formula": row["Formula"],
        "category": row["RateCategory"],
        "rate_multiplier": row["RateMultiplier"],
        "rate_sites": row["RateSites"],
        "valid_sites": row["ValidSites"],
        "skipped_sites": row["SkippedSites"],
        "satC": row["satC"],
        "satVar": row["satVar"],
        "satSE": row["satSE"],
        "satZ": row["satZ"],
        "weights": row["Weights"],
    }


def simulate_branch_batch(
    iqtree,
    work_root,
    model,
    sites,
    branch_length,
    branch_index,
    replicate_count,
    base_seed,
    keep_runs,
):
    """Generate all replicates for one branch in one independent SPRNG stream."""
    workdir = work_root / f"b{branch_tag(branch_length)}"
    workdir.mkdir(parents=True, exist_ok=True)
    tree = workdir / "tree.nwk"
    tree.write_text(f"((A:0.05,B:0.05):{branch_length:.12g},C:0.05,D:0.05);\n", encoding="utf-8")
    sim_prefix = workdir / "sim"
    seed = simulation_seed(base_seed, sites, branch_index, 0)
    command = [
        iqtree,
        "--alisim",
        str(sim_prefix),
        "-t",
        str(tree),
        "-m",
        model,
        "--length",
        str(sites),
        "--seed",
        str(seed),
        "--num-alignments",
        str(replicate_count),
        "-af",
        "fasta",
        "--skip-bl-check",
        "--redo",
        "--quiet",
    ]
    run_command(command)
    tasks = []
    for replicate in range(1, replicate_count + 1):
        alignment = (
            Path(str(sim_prefix) + ".fa")
            if replicate_count == 1
            else Path(str(sim_prefix) + f"_{replicate}.fa")
        )
        if not alignment.exists():
            raise FileNotFoundError(f"ALISIM batch output is missing: {alignment}")
        tasks.append(
            (
                iqtree,
                workdir,
                tree,
                alignment,
                model,
                branch_length,
                replicate,
                seed,
                keep_runs,
            )
        )
    return tasks


def analyze_batch_replicate(task):
    (
        iqtree,
        workdir,
        tree,
        alignment,
        model,
        branch_length,
        replicate,
        seed,
        keep_runs,
    ) = task
    sat_prefix = workdir / f"sat_r{replicate:04d}"
    run_command(
        [
            iqtree,
            "-s",
            str(alignment),
            "-te",
            str(tree),
            "-blfix",
            "-m",
            model,
            "--satute",
            "--prefix",
            str(sat_prefix),
            "-T",
            "1",
            "--redo",
            "--quiet",
        ]
    )
    pooled, all_target_rows = parse_stat(Path(str(sat_prefix) + ".sat.stat"))
    categories = parse_rate_categories(Path(str(sat_prefix) + ".iqtree"))
    rows = [pooled_to_detail(branch_length, replicate, seed, formula, pooled[formula]) for formula in FORMULAS]
    category_rows = [
        native_category_to_detail(branch_length, replicate, seed, row)
        for row in all_target_rows
        if row["RateCategory"] != "pooled" and row["Formula"] in {"dominant", "eigenvalue_weighted"}
    ]
    relative_weight_error = 0.0
    dominant_weight_error = 0.0
    for row in all_target_rows:
        if (
            row["RateCategory"] == "pooled"
            or row["Weights"] in {"", "NA"}
            or row["Weights"].startswith("soft_mixture_")
        ):
            continue
        weights = parse_numbers(row["Weights"])
        if row["Formula"] == "eigenvalue_weighted" and weights:
            relative_weight_error = max(relative_weight_error, abs(max(weights) - 1.0))
        if row["Formula"] == "dominant" and weights:
            dominant_weight_error = max(dominant_weight_error, *(abs(value - 1.0) for value in weights))
    metadata = {
        "categories": categories,
        "eigenvalues": parse_numbers(pooled["mixture_likelihood_weighted"]["Eigenvalues"]),
        "relative_weight_error": relative_weight_error,
        "dominant_weight_error": dominant_weight_error,
    }
    if not keep_runs:
        alignment.unlink(missing_ok=True)
        for output in workdir.glob(sat_prefix.name + "*"):
            if output.is_file():
                output.unlink()
    return rows, category_rows, metadata


def quantile(values, probability):
    ordered = sorted(values)
    if not ordered:
        return math.nan
    if len(ordered) == 1:
        return ordered[0]
    position = probability * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def mean_ci(values):
    mean = statistics.fmean(values)
    if len(values) < 2:
        return mean, mean, mean
    half = 1.96 * statistics.stdev(values) / math.sqrt(len(values))
    return mean, mean - half, mean + half


def wilson_interval(successes, total, z=1.96):
    if total <= 0:
        return math.nan, math.nan
    fraction = successes / total
    denominator = 1.0 + z * z / total
    center = (fraction + z * z / (2.0 * total)) / denominator
    half = z * math.sqrt(fraction * (1.0 - fraction) / total + z * z / (4.0 * total * total)) / denominator
    return center - half, center + half


def write_tsv(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def summarize(detail, null_length, null_reps):
    numeric = []
    for row in detail:
        converted = dict(row)
        for key in (
            "branch_length",
            "satC",
            "satVar",
            "satSE",
            "satZ",
            "satP",
            "information_fraction",
            "saturation_index",
        ):
            converted[key] = float(row[key])
        converted["replicate"] = int(row["replicate"])
        numeric.append(converted)

    training_limit = null_reps // 2
    thresholds = {}
    for formula in FORMULAS:
        training = [
            row["satZ"]
            for row in numeric
            if row["formula"] == formula
            and math.isclose(row["branch_length"], null_length)
            and row["replicate"] <= training_limit
        ]
        thresholds[formula] = quantile(training, 0.95)

    groups = defaultdict(list)
    for row in numeric:
        groups[(row["branch_length"], row["formula"])].append(row)

    summary = []
    for (branch_length, formula), rows in sorted(groups.items()):
        descriptive = rows
        decisions = rows
        if math.isclose(branch_length, null_length):
            decisions = [row for row in rows if row["replicate"] > training_limit]
        satc = [row["satC"] for row in descriptive]
        satse = [row["satSE"] for row in descriptive]
        satz = [row["satZ"] for row in descriptive]
        c_mean, c_low, c_high = mean_ci(satc)
        se_mean, se_low, se_high = mean_ci(satse)
        calibrated_success = sum(row["satZ"] > thresholds[formula] for row in decisions)
        nominal_success = sum(row["decision"] == "informative" for row in decisions)
        cal_low, cal_high = wilson_interval(calibrated_success, len(decisions))
        nom_low, nom_high = wilson_interval(nominal_success, len(decisions))
        summary.append(
            {
                "branch_length": f"{branch_length:.10g}",
                "formula": formula,
                "descriptive_replicates": len(descriptive),
                "decision_replicates": len(decisions),
                "satC_mean": c_mean,
                "satC_ci_low": c_low,
                "satC_ci_high": c_high,
                "satSE_mean": se_mean,
                "satSE_ci_low": se_low,
                "satSE_ci_high": se_high,
                "satZ_mean": statistics.fmean(satz),
                "satZ_median": quantile(satz, 0.5),
                "satZ_q10": quantile(satz, 0.1),
                "satZ_q90": quantile(satz, 0.9),
                "satZ_sd": statistics.stdev(satz) if len(satz) > 1 else math.nan,
                "information_fraction": statistics.fmean(
                    row["information_fraction"] for row in descriptive
                ),
                "saturation_index": statistics.fmean(
                    row["saturation_index"] for row in descriptive
                ),
                "nominal_fraction": nominal_success / len(decisions),
                "nominal_ci_low": nom_low,
                "nominal_ci_high": nom_high,
                "calibrated_threshold": thresholds[formula],
                "calibrated_fraction": calibrated_success / len(decisions),
                "calibrated_ci_low": cal_low,
                "calibrated_ci_high": cal_high,
            }
        )
    return summary, thresholds, numeric


def theoretical_rows(categories, eigenvalues, minimum, maximum, points=240):
    if minimum <= 0.0 or maximum <= minimum:
        raise ValueError("The theoretical curve requires 0 < minimum < maximum")
    nonzero = [value for value in eigenvalues if value < -1e-10]
    if not nonzero:
        raise ValueError("No negative nonstationary eigenvalues were found")
    dominant = max(nonzero)
    grid = [10 ** (math.log10(minimum) + index * (math.log10(maximum) - math.log10(minimum)) / (points - 1)) for index in range(points)]
    information_zero = sum(item["proportion"] for item in categories) * len(nonzero)
    theory = []
    for branch_length in grid:
        information = 0.0
        row = {"branch_length": branch_length}
        for item in categories:
            persistence = math.exp(dominant * item["rate"] * branch_length)
            row[f"rho_{item['category']}"] = persistence
            information += item["proportion"] * sum(
                math.exp(2.0 * eigenvalue * item["rate"] * branch_length) for eigenvalue in nonzero
            )
        row["information_fraction"] = information / information_zero
        row["saturation_index"] = 1.0 - row["information_fraction"]
        theory.append(row)
    return theory


def summarize_categories(category_detail, sites):
    numeric = []
    for row in category_detail:
        converted = dict(row)
        converted["branch_length"] = float(row["branch_length"])
        converted["replicate"] = int(row["replicate"])
        for key in ("rate_multiplier", "rate_sites", "valid_sites", "skipped_sites", "satC", "satVar", "satSE", "satZ"):
            converted[key] = float(row[key]) if row[key] not in {"", "NA"} else math.nan
        numeric.append(converted)
    groups = defaultdict(list)
    for row in numeric:
        groups[(row["branch_length"], row["formula"], row["category"])].append(row)
    summary = []
    for (branch_length, formula, category), rows in sorted(groups.items()):
        assignment = [row["rate_sites"] / sites for row in rows]
        valid_rows = [row for row in rows if math.isfinite(row["satC"]) and row["valid_sites"] > 0]
        contribution = [row["valid_sites"] / sites * row["satC"] for row in valid_rows]
        summary.append(
            {
                "branch_length": f"{branch_length:.10g}",
                "formula": formula,
                "category": category,
                "rate_multiplier": statistics.fmean(row["rate_multiplier"] for row in rows),
                "mean_assigned_fraction": statistics.fmean(assignment),
                "mean_category_satC": statistics.fmean(row["satC"] for row in valid_rows) if valid_rows else math.nan,
                "mean_pooled_contribution": statistics.fmean(contribution) if contribution else 0.0,
                "replicates": len(rows),
                "replicates_with_sites": len(valid_rows),
            }
        )
    return summary, numeric


def scale_linear(value, lower, upper, start, length, invert=False):
    if upper <= lower:
        return start + length / 2.0
    fraction = (value - lower) / (upper - lower)
    if invert:
        fraction = 1.0 - fraction
    return start + fraction * length


def scale_log_x(value, lower, upper, start, length):
    return scale_linear(math.log10(value), math.log10(lower), math.log10(upper), start, length)


def path(points):
    if not points:
        return ""
    return "M " + " L ".join(f"{x:.2f} {y:.2f}" for x, y in points)


def polygon(points):
    return " ".join(f"{x:.2f},{y:.2f}" for x, y in points)


def marker(svg, kind, x, y, color):
    if kind == "square":
        svg.append(f'<rect x="{x-3:.2f}" y="{y-3:.2f}" width="6" height="6" fill="white" stroke="{color}" stroke-width="1.7"/>')
    elif kind == "triangle":
        points = [(x, y - 3.8), (x - 3.5, y + 3.0), (x + 3.5, y + 3.0)]
        svg.append(f'<polygon points="{polygon(points)}" fill="white" stroke="{color}" stroke-width="1.7"/>')
    else:
        svg.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3" fill="white" stroke="{color}" stroke-width="1.7"/>')


def draw_axes(svg, left, top, width, height, title, subtitle, xmin, xmax, ymin, ymax, ylabel, log_y=False):
    svg.append(f'<text x="{left}" y="{top-28}" font-size="16" font-weight="700" fill="#202124">{escape(title)}</text>')
    svg.append(f'<text x="{left}" y="{top-10}" font-size="10.5" fill="#5f6368">{escape(subtitle)}</text>')
    svg.append(f'<rect x="{left}" y="{top}" width="{width}" height="{height}" fill="#ffffff" stroke="#3c4043" stroke-width="1"/>')
    for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
        y = top + (1.0 - fraction) * height
        raw = ymin + fraction * (ymax - ymin)
        label = f"10^{raw:.0f}" if log_y else f"{raw:.2g}"
        svg.append(f'<line x1="{left}" x2="{left+width}" y1="{y:.2f}" y2="{y:.2f}" stroke="#e8eaed"/>')
        svg.append(f'<text x="{left-7}" y="{y+4:.2f}" text-anchor="end" font-size="9.5" fill="#5f6368">{label}</text>')
    for tick in (0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0):
        if tick < xmin or tick > xmax:
            continue
        x = scale_log_x(tick, xmin, xmax, left, width)
        svg.append(f'<line x1="{x:.2f}" x2="{x:.2f}" y1="{top}" y2="{top+height}" stroke="#f1f3f4"/>')
        svg.append(f'<text x="{x:.2f}" y="{top+height+15}" text-anchor="middle" font-size="9.5" fill="#5f6368">{tick:g}</text>')
    svg.append(f'<text x="{left+width/2:.2f}" y="{top+height+34}" text-anchor="middle" font-size="10.5" fill="#3c4043">target branch length (log scale)</text>')
    x = left - 46
    y = top + height / 2
    svg.append(f'<text x="{x}" y="{y}" text-anchor="middle" font-size="10.5" fill="#3c4043" transform="rotate(-90 {x} {y})">{escape(ylabel)}</text>')


def render_figure(summary, theory, categories, sites, reps, null_reps, output_prefix):
    by_formula = defaultdict(list)
    for row in summary:
        converted = {key: (float(value) if key not in {"formula"} else value) for key, value in row.items() if key not in {"descriptive_replicates", "decision_replicates"}}
        by_formula[row["formula"]].append(converted)
    for rows in by_formula.values():
        rows.sort(key=lambda item: item["branch_length"])

    xmin = min(row["branch_length"] for rows in by_formula.values() for row in rows)
    xmax = max(row["branch_length"] for rows in by_formula.values() for row in rows)
    width, height = 1280, 1260
    panel_w, panel_h = 500, 275
    lefts = (105, 700)
    tops = (145, 515, 885)
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="640" y="38" text-anchor="middle" font-size="25" font-weight="700" fill="#202124">SatuTe movement toward saturation</text>',
        f'<text x="640" y="65" text-anchor="middle" font-size="12.5" fill="#5f6368">Fixed GTR+G4 model and tree; {sites:,} sites; {reps} alternative and {null_reps} null replicates</text>',
    ]

    plot_specs = [
        ("A  Raw coherence", "Mean across replicates; band is 95% CI", "satC_mean", "satC_ci_low", "satC_ci_high", "raw coherence, C̄", False),
        ("B  Standard error", "Mean across replicates; logarithmic y-axis", "satSE_mean", "satSE_ci_low", "satSE_ci_high", "satSE", True),
        ("C  Standardized statistic", "Median; band is the 10th–90th percentile", "satZ_median", "satZ_q10", "satZ_q90", "satZ", False),
        ("E  Calibrated informative probability", "Formula-specific null threshold; held-out null evaluation", "calibrated_fraction", "calibrated_ci_low", "calibrated_ci_high", "probability informative", False),
    ]

    positions = [(lefts[0], tops[0]), (lefts[1], tops[0]), (lefts[0], tops[1]), (lefts[0], tops[2])]
    for spec, (left, top) in zip(plot_specs, positions):
        title, subtitle, center_key, low_key, high_key, ylabel, log_y = spec
        all_values = []
        for rows in by_formula.values():
            for row in rows:
                for key in (low_key, high_key):
                    value = row[key]
                    if log_y:
                        if value > 0.0:
                            all_values.append(math.log10(value))
                    else:
                        all_values.append(value)
        if title.startswith("E"):
            ymin, ymax = 0.0, 1.0
        else:
            ymin, ymax = min(all_values), max(all_values)
            padding = 0.08 * max(ymax - ymin, 1e-9)
            ymin -= padding
            ymax += padding
            if not log_y and center_key == "satC_mean":
                ymin = min(ymin, 0.0)
            if not log_y and center_key == "satZ_median":
                ymin = min(ymin, 0.0)
        draw_axes(svg, left, top, panel_w, panel_h, title, subtitle, xmin, xmax, ymin, ymax, ylabel, log_y)
        if center_key == "satZ_median":
            y = scale_linear(1.6448536269514722, ymin, ymax, top, panel_h, invert=True)
            svg.append(f'<line x1="{left}" x2="{left+panel_w}" y1="{y:.2f}" y2="{y:.2f}" stroke="#3c4043" stroke-dasharray="5,4"/>')
            svg.append(f'<text x="{left+panel_w-4}" y="{y-5:.2f}" text-anchor="end" font-size="9" fill="#3c4043">nominal 5% threshold</text>')
        if title.startswith("E"):
            y = scale_linear(0.05, ymin, ymax, top, panel_h, invert=True)
            svg.append(f'<line x1="{left}" x2="{left+panel_w}" y1="{y:.2f}" y2="{y:.2f}" stroke="#3c4043" stroke-dasharray="5,4"/>')
        for formula in FORMULAS:
            label, color, dash, kind = FORMULA_STYLE[formula]
            rows = by_formula[formula]
            def transformed(value):
                if log_y:
                    return math.log10(max(value, 1e-300))
                return value
            upper = [(scale_log_x(row["branch_length"], xmin, xmax, left, panel_w), scale_linear(transformed(row[high_key]), ymin, ymax, top, panel_h, invert=True)) for row in rows]
            lower = [(scale_log_x(row["branch_length"], xmin, xmax, left, panel_w), scale_linear(transformed(row[low_key]), ymin, ymax, top, panel_h, invert=True)) for row in reversed(rows)]
            svg.append(f'<polygon points="{polygon(upper+lower)}" fill="{color}" fill-opacity="0.10" stroke="none"/>')
            points = [(scale_log_x(row["branch_length"], xmin, xmax, left, panel_w), scale_linear(transformed(row[center_key]), ymin, ymax, top, panel_h, invert=True)) for row in rows]
            dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
            svg.append(f'<path d="{path(points)}" fill="none" stroke="{color}" stroke-width="2.2"{dash_attr}/>')
            for x, y in points:
                marker(svg, kind, x, y, color)

    # D: category-specific persistence.
    left, top = lefts[1], tops[1]
    ymin, ymax = -12.0, 0.0
    draw_axes(svg, left, top, panel_w, panel_h, "D  Category-specific persistence", "Dominant mode; values below 10⁻¹² are clipped", xmin, xmax, ymin, ymax, "dominant persistence", True)
    category_colors = ("#0A3FBA", "#4778CE", "#829FE0", "#B7C7EF")
    for index, category in enumerate(categories):
        key = f"rho_{category['category']}"
        points = []
        for row in theory:
            value = max(row[key], 1e-12)
            points.append((scale_log_x(row["branch_length"], xmin, xmax, left, panel_w), scale_linear(math.log10(value), ymin, ymax, top, panel_h, invert=True)))
        color = category_colors[index % len(category_colors)]
        svg.append(f'<path d="{path(points)}" fill="none" stroke="{color}" stroke-width="2.1" stroke-dasharray="{2+index*3},{3+index}"/>')
        legend_x = left + 12
        legend_y = top + 17 + index * 15
        svg.append(f'<line x1="{legend_x}" x2="{legend_x+25}" y1="{legend_y}" y2="{legend_y}" stroke="{color}" stroke-width="2.1" stroke-dasharray="{2+index*3},{3+index}"/>')
        svg.append(f'<text x="{legend_x+31}" y="{legend_y+3.5}" font-size="9.5" fill="{color}">c{category["category"]}: r={category["rate"]:.3g}</text>')

    # F: monotonic model-based saturation index.
    left, top = lefts[1], tops[2]
    draw_axes(svg, left, top, panel_w, panel_h, "F  Model-based saturation index", "Equal spectral coefficients, a_k=1; independent of alignment size", xmin, xmax, 0.0, 1.0, "1 − I(t)/I(0)", False)
    points = [(scale_log_x(row["branch_length"], xmin, xmax, left, panel_w), scale_linear(row["saturation_index"], 0.0, 1.0, top, panel_h, invert=True)) for row in theory]
    svg.append(f'<path d="{path(points)}" fill="none" stroke="#202124" stroke-width="2.6"/>')
    native_points = [
        (
            scale_log_x(row["branch_length"], xmin, xmax, left, panel_w),
            scale_linear(row["saturation_index"], 0.0, 1.0, top, panel_h, invert=True),
        )
        for row in by_formula["dominant"]
    ]
    for x, y in native_points:
        svg.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3.2" fill="#ffffff" stroke="#202124" stroke-width="1.5"/>')

    # Shared legend.
    legend_y = 98
    for index, formula in enumerate(FORMULAS):
        label, color, dash, kind = FORMULA_STYLE[formula]
        x = 175 + index * 330
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        svg.append(f'<line x1="{x}" x2="{x+38}" y1="{legend_y}" y2="{legend_y}" stroke="{color}" stroke-width="2.5"{dash_attr}/>')
        marker(svg, kind, x + 19, legend_y, color)
        svg.append(f'<text x="{x+48}" y="{legend_y+4}" font-size="11" fill="#3c4043">{escape(label)}</text>')

    svg.append('<text x="105" y="1240" font-size="10.5" fill="#5f6368">Raw coherence is formula-scale dependent. satZ is a detection statistic; the saturation index is the monotonic model-based ruler.</text>')
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
    return svg_path, pdf_path if pdf_path.exists() else None, png_path if png_path.exists() else None


def render_category_figure(category_summary, summary, output_prefix):
    rows = []
    for row in category_summary:
        if row["formula"] != "eigenvalue_weighted":
            continue
        converted = dict(row)
        for key in ("branch_length", "rate_multiplier", "mean_assigned_fraction", "mean_category_satC", "mean_pooled_contribution"):
            converted[key] = float(row[key])
        rows.append(converted)
    categories = sorted({row["category"] for row in rows}, key=int)
    by_category = defaultdict(list)
    for row in rows:
        by_category[row["category"]].append(row)
    for values in by_category.values():
        values.sort(key=lambda row: row["branch_length"])
    pooled = sorted(
        (
            (float(row["branch_length"]), float(row["satC_mean"]))
            for row in summary
            if row["formula"] == "eigenvalue_weighted"
        ),
        key=lambda item: item[0],
    )
    xmin = min(row["branch_length"] for row in rows)
    xmax = max(row["branch_length"] for row in rows)
    contribution_xmin = 10.0
    width, height = 1180, 520
    panel_w, panel_h = 450, 300
    lefts = (95, 670)
    top = 140
    category_colors = ("#0A3FBA", "#4778CE", "#829FE0", "#B7C7EF")
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="590" y="36" text-anchor="middle" font-size="24" font-weight="700" fill="#202124">Why relative-weighted coherence changes abruptly</text>',
        '<text x="590" y="62" text-anchor="middle" font-size="12" fill="#5f6368">Hard MAP rate categories; fixed GTR+G4 model and 1,000-site alignments</text>',
    ]
    draw_axes(svg, lefts[0], top, panel_w, panel_h, "A  Hard-category occupancy", "Mean fraction of sites assigned to each category", xmin, xmax, 0.0, 1.0, "assigned fraction", False)
    contribution_values = [
        row["mean_pooled_contribution"] for row in rows if row["branch_length"] >= contribution_xmin
    ] + [value for branch, value in pooled if branch >= contribution_xmin]
    ymin, ymax = min(contribution_values), max(contribution_values)
    padding = 0.08 * (ymax - ymin)
    ymin -= padding
    ymax += padding
    draw_axes(svg, lefts[1], top, panel_w, panel_h, "B  Signed contribution to pooled coherence", "Transition region t≥10; contribution = (valid sites / all sites) × category satC", contribution_xmin, xmax, ymin, ymax, "coherence contribution", False)
    zero_y = scale_linear(0.0, ymin, ymax, top, panel_h, invert=True)
    svg.append(f'<line x1="{lefts[1]}" x2="{lefts[1]+panel_w}" y1="{zero_y:.2f}" y2="{zero_y:.2f}" stroke="#3c4043" stroke-dasharray="5,4"/>')

    for index, category in enumerate(categories):
        color = category_colors[index % len(category_colors)]
        values = by_category[category]
        dash = f"{2+index*3},{3+index}"
        occupancy_points = [
            (
                scale_log_x(row["branch_length"], xmin, xmax, lefts[0], panel_w),
                scale_linear(row["mean_assigned_fraction"], 0.0, 1.0, top, panel_h, invert=True),
            )
            for row in values
        ]
        contribution_points = [
            (
                scale_log_x(row["branch_length"], contribution_xmin, xmax, lefts[1], panel_w),
                scale_linear(row["mean_pooled_contribution"], ymin, ymax, top, panel_h, invert=True),
            )
            for row in values
            if row["branch_length"] >= contribution_xmin
        ]
        svg.append(f'<path d="{path(occupancy_points)}" fill="none" stroke="{color}" stroke-width="2.2" stroke-dasharray="{dash}"/>')
        svg.append(f'<path d="{path(contribution_points)}" fill="none" stroke="{color}" stroke-width="2.2" stroke-dasharray="{dash}"/>')
        legend_x = 170 + index * 225
        svg.append(f'<line x1="{legend_x}" x2="{legend_x+34}" y1="84" y2="84" stroke="{color}" stroke-width="2.3" stroke-dasharray="{dash}"/>')
        rate = values[0]["rate_multiplier"]
        svg.append(f'<text x="{legend_x+42}" y="88" font-size="10.5" fill="{color}">category {category}, r={rate:.3g}</text>')

    pooled_points = [
        (scale_log_x(branch, contribution_xmin, xmax, lefts[1], panel_w), scale_linear(value, ymin, ymax, top, panel_h, invert=True))
        for branch, value in pooled
        if branch >= contribution_xmin
    ]
    svg.append(f'<path d="{path(pooled_points)}" fill="none" stroke="#202124" stroke-width="2.7"/>')
    svg.append(f'<text x="{lefts[1]+panel_w-6}" y="{top+18}" text-anchor="end" font-size="10" fill="#202124">black: pooled total</text>')
    svg.append('<text x="95" y="505" font-size="10.5" fill="#5f6368">Abrupt pooled movement is caused by reassignment and cancellation between signed category contributions, not by discontinuous eigenvalue weights.</text>')
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
    return svg_path, pdf_path if pdf_path.exists() else None, png_path if png_path.exists() else None


def check_report(
    detail,
    category_detail,
    summary,
    thresholds,
    categories,
    eigenvalues,
    theory,
    null_length,
    normalization_checks,
    output,
):
    checks = []
    z_errors = []
    for row in detail:
        satc, satse, satz = float(row["satC"]), float(row["satSE"]), float(row["satZ"])
        if satse > 0.0:
            z_errors.append(abs(satc / satse - satz))
    checks.append(("Arithmetic identity", "PASS" if max(z_errors, default=0.0) < 1e-7 else "FAIL", f"max |satC/satSE − satZ| = {max(z_errors, default=0.0):.3g}"))

    dominant_error = normalization_checks.get("dominant_weights_equal_one", math.inf)
    relative_error = normalization_checks.get("relative_maximum_weight_equals_one", math.inf)
    checks.append(("Dominant-weight normalization", "PASS" if dominant_error < 1e-12 else "FAIL", f"max |w − 1| = {dominant_error:.3g}"))
    checks.append(("Relative-weight normalization", "PASS" if relative_error < 1e-12 else "FAIL", f"max |max_k w_k − 1| = {relative_error:.3g}"))

    prop_sum = sum(item["proportion"] for item in categories)
    rates_valid = all(item["rate"] >= 0.0 and item["proportion"] > 0.0 for item in categories)
    checks.append(("Rate-category contract", "PASS" if rates_valid and abs(prop_sum - 1.0) < 1e-10 else "FAIL", f"sum q_c = {prop_sum:.12g}; rates nonnegative = {rates_valid}"))

    monotonic_persistence = True
    for category in categories:
        values = [row[f"rho_{category['category']}"] for row in theory]
        monotonic_persistence &= all(right <= left + 1e-14 for left, right in zip(values, values[1:]))
    checks.append(("Persistence monotonicity", "PASS" if monotonic_persistence else "FAIL", "every category-specific dominant persistence curve is nonincreasing"))

    saturation = [row["saturation_index"] for row in theory]
    monotonic_saturation = all(right + 1e-14 >= left for left, right in zip(saturation, saturation[1:]))
    checks.append(("Saturation-index monotonicity", "PASS" if monotonic_saturation else "FAIL", f"S(0)=0 analytically; plotted S range = {saturation[0]:.6g} to {saturation[-1]:.6g}"))

    nonzero_eigenvalues = [value for value in eigenvalues if value < -1e-10]
    native_scale_errors = []
    native_scale_by_branch = defaultdict(list)
    for row in detail:
        branch_length = float(row["branch_length"])
        expected_information = sum(
            category["proportion"]
            * sum(
                math.exp(2.0 * eigenvalue * category["rate"] * branch_length)
                for eigenvalue in nonzero_eigenvalues
            )
            / len(nonzero_eigenvalues)
            for category in categories
        ) / sum(category["proportion"] for category in categories)
        observed_information = float(row["information_fraction"])
        observed_saturation = float(row["saturation_index"])
        native_scale_errors.extend(
            (
                abs(observed_information - expected_information),
                abs(observed_saturation - (1.0 - expected_information)),
                abs(observed_information + observed_saturation - 1.0),
            )
        )
        native_scale_by_branch[branch_length].append(observed_saturation)
    formula_invariance_error = max(
        (max(values) - min(values) for values in native_scale_by_branch.values()),
        default=math.inf,
    )
    max_native_scale_error = max(native_scale_errors, default=math.inf)
    # The independently parsed .iqtree table prints category rates to four
    # significant digits, whereas native SatuTe uses the full-precision rates.
    scale_status = "PASS" if max_native_scale_error < 2e-5 and formula_invariance_error < 1e-12 else "FAIL"
    checks.append(("Native saturation-scale identity", scale_status, f"max equation error = {max_native_scale_error:.3g} with rounded report rates; max between-formula difference = {formula_invariance_error:.3g}"))

    pooled_by_key = {
        (float(row["branch_length"]), int(row["replicate"]), row["formula"]): row
        for row in detail
        if row["formula"] in {"dominant", "eigenvalue_weighted"}
    }
    categories_by_key = defaultdict(list)
    for row in category_detail:
        categories_by_key[(float(row["branch_length"]), int(row["replicate"]), row["formula"])].append(row)
    pooled_mean_errors = []
    pooled_variance_errors = []
    for key, rows in categories_by_key.items():
        valid = [row for row in rows if row["valid_sites"] not in {"", "0"} and row["satC"] not in {"", "NA"}]
        total = sum(float(row["valid_sites"]) for row in valid)
        if total <= 0.0 or key not in pooled_by_key:
            continue
        recomputed_mean = sum(float(row["valid_sites"]) / total * float(row["satC"]) for row in valid)
        recomputed_variance = sum(float(row["valid_sites"]) / total * float(row["satVar"]) for row in valid)
        pooled_mean_errors.append(abs(recomputed_mean - float(pooled_by_key[key]["satC"])))
        pooled_variance_errors.append(abs(recomputed_variance - float(pooled_by_key[key]["satVar"])))
    max_pool_mean_error = max(pooled_mean_errors, default=math.inf)
    max_pool_variance_error = max(pooled_variance_errors, default=math.inf)
    pool_status = "PASS" if max_pool_mean_error < 1e-7 and max_pool_variance_error < 1e-7 else "FAIL"
    checks.append(("Hard-category pooling identity", pool_status, f"max mean error = {max_pool_mean_error:.3g}; max variance error = {max_pool_variance_error:.3g}"))

    assignment_by_length = defaultdict(lambda: defaultdict(list))
    for row in category_detail:
        if row["formula"] != "dominant":
            continue
        assignment_by_length[float(row["branch_length"])][row["category"]].append(float(row["rate_sites"]))
    assignment_composition = {
        branch_length: {
            category: statistics.fmean(counts) / sum(statistics.fmean(values) for values in categories_at_length.values())
            for category, counts in categories_at_length.items()
        }
        for branch_length, categories_at_length in assignment_by_length.items()
    }
    largest_shift = (0.0, math.nan, math.nan)
    ordered_lengths = sorted(assignment_composition)
    for left_length, right_length in zip(ordered_lengths, ordered_lengths[1:]):
        category_names = set(assignment_composition[left_length]) | set(assignment_composition[right_length])
        total_variation = 0.5 * sum(
            abs(assignment_composition[left_length].get(category, 0.0) - assignment_composition[right_length].get(category, 0.0))
            for category in category_names
        )
        if total_variation > largest_shift[0]:
            largest_shift = (total_variation, left_length, right_length)
    checks.append(("Largest hard-assignment shift", "OBSERVED", f"total-variation change = {largest_shift[0]:.3f} between t={largest_shift[1]:g} and t={largest_shift[2]:g}"))

    summary_by_key = {(float(row["branch_length"]), row["formula"]): row for row in summary}
    branch_lengths = sorted({float(row["branch_length"]) for row in summary})
    first = branch_lengths[0]
    for formula in FORMULAS:
        start = float(summary_by_key[(first, formula)]["satC_mean"])
        end = float(summary_by_key[(null_length, formula)]["satC_mean"])
        status = "PASS" if abs(end) < abs(start) else "FAIL"
        checks.append((f"Raw-coherence endpoint: {formula}", status, f"|mean satC|: {abs(start):.4g} → {abs(end):.4g}"))

    dominant_end = float(summary_by_key[(null_length, "dominant")]["satC_mean"])
    weighted_end = float(summary_by_key[(null_length, "eigenvalue_weighted")]["satC_mean"])
    convergence_error = abs(dominant_end - weighted_end)
    checks.append(("Long-branch weighted-to-dominant convergence", "PASS" if convergence_error < 1e-8 else "FAIL", f"difference in mean satC at t={null_length:g}: {convergence_error:.3g}"))

    null_rows = [row for row in summary if math.isclose(float(row["branch_length"]), null_length)]
    for row in null_rows:
        mean_z = float(row["satZ_mean"])
        sd_z = float(row["satZ_sd"])
        status = "PASS" if abs(mean_z) <= 0.2 and 0.8 <= sd_z <= 1.2 else "FAIL"
        checks.append((f"Nominal normal null: {row['formula']}", status, f"mean Z = {mean_z:.3f}; SD = {sd_z:.3f}; q95 training = {thresholds[row['formula']]:.3f}"))
        fraction = float(row["calibrated_fraction"])
        low = float(row["calibrated_ci_low"])
        high = float(row["calibrated_ci_high"])
        status = "PASS" if low <= 0.05 <= high else "FAIL"
        checks.append((f"Held-out calibrated null: {row['formula']}", status, f"fraction = {fraction:.3f}; Wilson 95% CI {low:.3f}–{high:.3f}"))

    lines = [
        "# Saturation diagnostic assumption checks",
        "",
        "These checks distinguish algebraic requirements from empirical observations. A failed normal-null check does not invalidate the statistic; it invalidates use of the uncalibrated standard-normal threshold in that scenario.",
        "",
        "| Step | Check | Status | Evidence |",
        "|---:|---|---|---|",
    ]
    for index, (name, status, evidence) in enumerate(checks, 1):
        safe_evidence = evidence.replace("|", "\\|")
        lines.append(f"| {index} | {name} | **{status}** | {safe_evidence} |")
    lines.extend(
        [
            "",
            "## Tested definition",
            "",
            "With equal nonstationary spectral coefficients,",
            "",
            "$$",
            "I(t)=\\sum_c q_c\\sum_{k\\ne0}\\exp(2\\lambda_k r_c t),",
            "\\qquad",
            "S(t)=1-\\frac{I(t)}{I(0)}.",
            "$$",
            "",
            "For $q_c>0$, $r_c\\ge0$, and $\\lambda_k<0$, every derivative term is nonpositive, so $I(t)$ is nonincreasing and $S(t)$ is nondecreasing.",
        ]
    )
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def load_detail(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_categories(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [{"category": row["category"], "rate": float(row["rate"]), "proportion": float(row["proportion"])} for row in rows]


def load_normalization_checks(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return {row["check"]: float(row["max_absolute_error"]) for row in csv.DictReader(handle, delimiter="\t")}


def main():
    parser = argparse.ArgumentParser(description="Simulate and plot native SatuTe saturation diagnostics.")
    parser.add_argument("--iqtree", default=str(IQTREE_DEFAULT))
    parser.add_argument(
        "--outdir", default=str(PROJECT_ROOT / "artifacts/local/saturation-diagnostics")
    )
    parser.add_argument(
        "--figure-prefix",
        default=str(PROJECT_ROOT / "artifacts/local/saturation-diagnostics/figures/main"),
    )
    parser.add_argument(
        "--category-figure-prefix",
        default=str(PROJECT_ROOT / "artifacts/local/saturation-diagnostics/figures/categories"),
    )
    parser.add_argument("--model", default=DEFAULT_MODEL)
    parser.add_argument("--sites", type=int, default=1000)
    parser.add_argument("--branch-lengths", default=DEFAULT_BRANCHES)
    parser.add_argument("--null-length", type=float, default=10000.0)
    parser.add_argument("--reps", type=int, default=100)
    parser.add_argument("--null-reps", type=int, default=500)
    parser.add_argument("--workers", type=int, default=max(1, min(8, os.cpu_count() or 1)))
    parser.add_argument("--base-seed", type=int, default=1200000)
    parser.add_argument("--keep-runs", action="store_true")
    parser.add_argument("--plot-only", action="store_true")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    detail_path = outdir / "saturation_diagnostic_detail.tsv"
    category_path = outdir / "rate_categories.tsv"
    eigenvalue_path = outdir / "eigenvalues.tsv"
    category_detail_path = outdir / "native_category_detail.tsv"

    if args.plot_only:
        detail = load_detail(detail_path)
        category_detail = load_detail(category_detail_path)
        categories = load_categories(category_path)
        normalization_checks = load_normalization_checks(outdir / "weight_normalization_checks.tsv")
        with eigenvalue_path.open("r", encoding="utf-8", newline="") as handle:
            eigenvalues = [float(row["eigenvalue"]) for row in csv.DictReader(handle, delimiter="\t")]
    else:
        iqtree = str(Path(args.iqtree).resolve())
        if not os.path.isfile(iqtree) or not os.access(iqtree, os.X_OK):
            raise SystemExit(f"Cannot execute IQ-TREE: {iqtree}")
        branch_lengths = sorted(set(parse_numbers(args.branch_lengths)))
        if args.null_length not in branch_lengths:
            branch_lengths.append(args.null_length)
            branch_lengths.sort()
        shutil.rmtree(outdir / "runs", ignore_errors=True)
        tasks = []
        for branch_index, branch_length in enumerate(branch_lengths):
            replicate_count = args.null_reps if math.isclose(branch_length, args.null_length) else args.reps
            tasks.extend(
                simulate_branch_batch(
                    iqtree,
                    outdir / "runs",
                    args.model,
                    args.sites,
                    branch_length,
                    branch_index,
                    replicate_count,
                    args.base_seed,
                    args.keep_runs,
                )
            )
        detail = []
        category_detail = []
        categories = None
        eigenvalues = None
        relative_weight_error = 0.0
        dominant_weight_error = 0.0
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = [executor.submit(analyze_batch_replicate, task) for task in tasks]
            completed = 0
            for future in as_completed(futures):
                rows, native_category_rows, metadata = future.result()
                detail.extend(rows)
                category_detail.extend(native_category_rows)
                completed += 1
                if completed % 100 == 0 or completed == len(tasks):
                    print(f"completed {completed}/{len(tasks)} replicate tasks", flush=True)
                if categories is None:
                    categories = metadata["categories"]
                    eigenvalues = metadata["eigenvalues"]
                relative_weight_error = max(relative_weight_error, metadata["relative_weight_error"])
                dominant_weight_error = max(dominant_weight_error, metadata["dominant_weight_error"])
        if not args.keep_runs:
            shutil.rmtree(outdir / "runs", ignore_errors=True)
        detail.sort(key=lambda row: (float(row["branch_length"]), int(row["replicate"]), row["formula"]))
        category_detail.sort(key=lambda row: (float(row["branch_length"]), int(row["replicate"]), row["formula"], int(row["category"])))
        write_tsv(detail_path, DETAIL_FIELDS, detail)
        write_tsv(category_detail_path, CATEGORY_FIELDS, category_detail)
        write_tsv(category_path, ("category", "rate", "proportion"), categories)
        write_tsv(eigenvalue_path, ("mode", "eigenvalue"), [{"mode": index, "eigenvalue": value} for index, value in enumerate(eigenvalues)])
        (outdir / "weight_normalization_checks.tsv").write_text(
            "check\tmax_absolute_error\n"
            f"dominant_weights_equal_one\t{dominant_weight_error:.12g}\n"
            f"relative_maximum_weight_equals_one\t{relative_weight_error:.12g}\n",
            encoding="utf-8",
        )
        normalization_checks = {
            "dominant_weights_equal_one": dominant_weight_error,
            "relative_maximum_weight_equals_one": relative_weight_error,
        }

    summary, thresholds, _numeric = summarize(detail, args.null_length, args.null_reps)
    summary_fields = tuple(summary[0].keys())
    write_tsv(outdir / "saturation_diagnostic_summary.tsv", summary_fields, summary)
    write_tsv(outdir / "calibrated_thresholds.tsv", ("formula", "q95"), [{"formula": formula, "q95": thresholds[formula]} for formula in FORMULAS])
    category_summary, _category_numeric = summarize_categories(category_detail, args.sites)
    write_tsv(outdir / "native_category_summary.tsv", tuple(category_summary[0].keys()), category_summary)
    minimum = min(float(row["branch_length"]) for row in detail)
    maximum = max(float(row["branch_length"]) for row in detail)
    theory = theoretical_rows(categories, eigenvalues, minimum, maximum)
    theory_fields = tuple(theory[0].keys())
    write_tsv(outdir / "model_saturation_curve.tsv", theory_fields, theory)
    figure_paths = render_figure(summary, theory, categories, args.sites, args.reps, args.null_reps, Path(args.figure_prefix))
    category_figure_paths = render_category_figure(category_summary, summary, Path(args.category_figure_prefix))
    check_report(
        detail,
        category_detail,
        summary,
        thresholds,
        categories,
        eigenvalues,
        theory,
        args.null_length,
        normalization_checks,
        outdir / "assumption_checks.md",
    )
    print(detail_path)
    print(outdir / "saturation_diagnostic_summary.tsv")
    for figure_path in figure_paths:
        if figure_path:
            print(figure_path)
    for figure_path in category_figure_paths:
        if figure_path:
            print(figure_path)


if __name__ == "__main__":
    main()
