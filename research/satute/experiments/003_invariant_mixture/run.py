#!/usr/bin/env python3
"""Paired exact-pattern simulation for the refined weighted +I statistic.

The four-taxon pattern distribution is enumerated exactly under fixed GTR+I+G4
and GTR+I+R4 models.  Multinomial alignments drawn from the same counts are
then evaluated with both the former counterfactual construction and the
refined construction.  This pairing isolates the formula change.
"""

import argparse
import csv
import hashlib
import itertools
import math
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from satute_analysis.models import build_q, reversible_eigendecomposition, transition_matrix
from satute_analysis.statistics import mean_interval, paired_bootstrap_interval, wilson_interval


SUBSTITUTION_MODEL = "GTR{1,2,1,1,2,1}+F{0.30,0.20,0.20,0.30}"
GAMMA_RATES = np.array([0.03338775338, 0.2519159176, 0.8202684819, 2.894427847])
FREE_PROPORTIONS = np.array([0.35, 0.25, 0.23, 0.17])
FREE_RATES = np.array([0.07476557822, 0.4208486335, 1.158502766, 3.542142665])
INVARIANT_PROPORTION = 0.10
NOMINAL_THRESHOLD = 1.6448536269514722
PATTERNS = tuple(itertools.product(range(4), repeat=4))
METHODS = ("legacy_counterfactual", "refined_invariant_null")
COLORS = {"legacy_counterfactual": "#C23B33", "refined_invariant_null": "#0A6F5A"}
DEFAULT_BRANCH_LENGTHS = (
    "0.5,1,2,3,4,5,7.5,10,12.5,15,17.5,20,22.5,25,30,35,40,45,50,"
    "60,75,100,150,200,500,1000,10000"
)


def model_grid():
    variable_mass = 1.0 - INVARIANT_PROPORTION
    return {
        "+I+G4": {
            "rates": np.concatenate(([0.0], GAMMA_RATES / variable_mass)),
            "proportions": np.concatenate(
                ([INVARIANT_PROPORTION], np.repeat(variable_mass / 4.0, 4))
            ),
            "iqtree_model": SUBSTITUTION_MODEL + "+I{0.1}+G4{0.5}",
        },
        "+I+R4": {
            "rates": np.concatenate(([0.0], FREE_RATES / variable_mass)),
            "proportions": np.concatenate(
                ([INVARIANT_PROPORTION], variable_mass * FREE_PROPORTIONS)
            ),
            "iqtree_model": SUBSTITUTION_MODEL
            + "+I{0.1}+R4{0.35,0.07476557822,0.25,0.4208486335,"
            "0.23,1.158502766,0.17,3.542142665}",
        },
    }


def category_terms(eigenvalues, eigenvectors, inverse, pi, modes, rate):
    transition = transition_matrix(eigenvalues, eigenvectors, inverse, 0.05 * rate)
    lefts, rights, independent, joint, coherence = [], [], [], [], []
    for pattern in PATTERNS:
        left = transition[:, pattern[0]] * transition[:, pattern[1]]
        right = transition[:, pattern[2]] * transition[:, pattern[3]]
        left_probability = float(pi @ left)
        right_probability = float(pi @ right)
        lefts.append(left)
        rights.append(right)
        independent.append(left_probability * right_probability)
        joint.append(float(pi @ (left * right)))
        if left_probability > 0.0 and right_probability > 0.0:
            left_posterior = pi * left / left_probability
            right_posterior = pi * right / right_probability
            coherence.append(
                (eigenvectors[:, modes].T @ left_posterior)
                * (eigenvectors[:, modes].T @ right_posterior)
            )
        else:
            coherence.append(np.zeros(len(modes)))
    return {
        "left": np.asarray(lefts),
        "right": np.asarray(rights),
        "independent": np.asarray(independent),
        "joint": np.asarray(joint),
        "coherence": np.asarray(coherence),
    }


def exact_cell(eigenvalues, eigenvectors, inverse, pi, modes, model, branch_length):
    rates = model["rates"]
    proportions = model["proportions"]
    categories = [
        category_terms(eigenvalues, eigenvectors, inverse, pi, modes, float(rate))
        for rate in rates
    ]
    refined_null = np.zeros(len(PATTERNS))
    counterfactual_null = np.zeros(len(PATTERNS))
    alternative = np.zeros(len(PATTERNS))
    refined_numerator = np.zeros(len(PATTERNS))
    counterfactual_numerator = np.zeros(len(PATTERNS))
    positive_exponents = [
        eigenvalues[mode] * branch_length * rate
        for rate, proportion in zip(rates, proportions)
        if rate > 0.0 and proportion > 0.0
        for mode in modes
    ]
    refined_shift = max(positive_exponents)

    for rate, proportion, category in zip(rates, proportions, categories):
        counterfactual_null += proportion * category["independent"]
        counterfactual_weights = np.exp(eigenvalues[modes] * branch_length * rate)
        counterfactual_numerator += (
            proportion
            * category["independent"]
            * (category["coherence"] @ counterfactual_weights)
        )

        if rate == 0.0:
            refined_null += proportion * category["joint"]
        else:
            refined_null += proportion * category["independent"]
            refined_weights = np.exp(
                eigenvalues[modes] * branch_length * rate - refined_shift
            )
            refined_numerator += (
                proportion
                * category["independent"]
                * (category["coherence"] @ refined_weights)
            )

        focal = transition_matrix(
            eigenvalues, eigenvectors, inverse, branch_length * rate
        )
        alternative += proportion * np.asarray(
            [
                np.einsum("a,a,ab,b->", pi, left, focal, right)
                for left, right in zip(category["left"], category["right"])
            ]
        )

    refined_null /= refined_null.sum()
    counterfactual_null /= counterfactual_null.sum()
    alternative /= alternative.sum()
    scores = {
        "legacy_counterfactual": counterfactual_numerator / counterfactual_null,
        "refined_invariant_null": refined_numerator / refined_null,
    }
    likelihood_ratio = alternative / refined_null
    identity_error = float(
        np.max(
            np.abs(
                likelihood_ratio
                - 1.0
                - math.exp(refined_shift) * scores["refined_invariant_null"]
            )
        )
    )
    return refined_null, alternative, scores, refined_shift, identity_error


def seed_for(base_seed, *parts):
    payload = ":".join(str(part) for part in (base_seed,) + parts).encode("ascii")
    value = int.from_bytes(hashlib.blake2s(payload, digest_size=8).digest(), "big")
    return 1 + value % 2147483646


def simulate_counts(probabilities, nsites, replicates, rng):
    return rng.multinomial(nsites, probabilities, size=replicates)


def studentized(counts, score, nsites):
    means = counts @ score / nsites
    centered_sum = counts @ (score * score) - nsites * means * means
    variances = np.maximum(centered_sum / (nsites - 1.0), 0.0)
    denominator = np.sqrt(variances / nsites)
    return np.divide(
        means,
        denominator,
        out=np.full_like(means, np.nan),
        where=denominator > 0.0,
    )


def write_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def render_figure(summary_rows, paired_rows, output_prefix, png_scale=3.0):
    summary = []
    for row in summary_rows:
        summary.append(
            {
                key: (value if key in {"model", "method"} else float(value))
                for key, value in row.items()
            }
        )
    paired = []
    for row in paired_rows:
        paired.append(
            {
                key: (value if key in {"model", "rule", "sample"} else float(value))
                for key, value in row.items()
            }
        )

    width, height = 1240, 850
    panel_w, panel_h = 480, 250
    panels = [(90, 155), (675, 155), (90, 545), (675, 545)]
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="620" y="38" text-anchor="middle" font-size="24" font-weight="700" fill="#202124">Weighted +I correction: paired exact-pattern simulation</text>',
        '<text x="620" y="64" text-anchor="middle" font-size="12" fill="#5f6368">Same multinomial alignment counts evaluated by both formulas; one-sided alpha = 0.05</text>',
    ]

    def axes(left, top, title, subtitle, xmin, xmax, ymin, ymax, xlabel, ylabel, log_x=False):
        svg.append(f'<text x="{left}" y="{top-28}" font-size="15" font-weight="700" fill="#202124">{title}</text>')
        svg.append(f'<text x="{left}" y="{top-10}" font-size="10.5" fill="#5f6368">{subtitle}</text>')
        svg.append(f'<rect x="{left}" y="{top}" width="{panel_w}" height="{panel_h}" fill="#fff" stroke="#3c4043"/>')
        for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
            y = top + panel_h * (1.0 - fraction)
            value = ymin + fraction * (ymax - ymin)
            svg.append(f'<line x1="{left}" x2="{left+panel_w}" y1="{y:.2f}" y2="{y:.2f}" stroke="#eceff1"/>')
            svg.append(f'<text x="{left-7}" y="{y+4:.2f}" text-anchor="end" font-size="9" fill="#5f6368">{value:.2g}</text>')
        ticks = (100, 500, 1000, 5000, 10000) if log_x else (0.5, 2, 5, 10, 20, 50, 200, 10000)
        for value in ticks:
            if value < xmin or value > xmax:
                continue
            fraction = ((math.log10(value)-math.log10(xmin))/(math.log10(xmax)-math.log10(xmin))) if log_x else ((math.log10(value)-math.log10(xmin))/(math.log10(xmax)-math.log10(xmin)))
            x = left + panel_w * fraction
            svg.append(f'<text x="{x:.2f}" y="{top+panel_h+16}" text-anchor="middle" font-size="8.5" fill="#5f6368">{value:g}</text>')
        svg.append(f'<text x="{left+panel_w/2}" y="{top+panel_h+36}" text-anchor="middle" font-size="10" fill="#3c4043">{xlabel}</text>')
        svg.append(f'<text x="{left-53}" y="{top+panel_h/2}" text-anchor="middle" font-size="10" fill="#3c4043" transform="rotate(-90 {left-53} {top+panel_h/2})">{ylabel}</text>')

    def polyline(rows, xkey, ykey, left, top, xmin, xmax, ymin, ymax, color, dash=""):
        points = []
        for row in sorted(rows, key=lambda item: item[xkey]):
            x = left + panel_w * (math.log10(row[xkey])-math.log10(xmin)) / (math.log10(xmax)-math.log10(xmin))
            y = top + panel_h * (1.0 - (row[ykey]-ymin)/(ymax-ymin))
            points.append((x, y))
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        svg.append('<polyline points="' + " ".join(f"{x:.2f},{y:.2f}" for x, y in points) + f'" fill="none" stroke="{color}" stroke-width="2.2"{dash_attr}/>')
        for x, y in points:
            svg.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2.7" fill="#fff" stroke="{color}" stroke-width="1.5"/>')

    # Extreme-branch nominal false positives by alignment length.
    left, top = panels[0]
    axes(left, top, "A  Extreme-branch false positives", "t = 10,000 under the correct +I null", 100, 10000, 0, 1, "alignment length, n (log scale)", "nominal rejection", True)
    y05 = top + panel_h * 0.95
    svg.append(f'<line x1="{left}" x2="{left+panel_w}" y1="{y05:.2f}" y2="{y05:.2f}" stroke="#3c4043" stroke-dasharray="4,4"/>')
    for model_index, model in enumerate(("+I+G4", "+I+R4")):
        for method in METHODS:
            rows = [row for row in summary if row["model"] == model and row["method"] == method and row["branch_length"] == 10000]
            polyline(rows, "nsites", "null_nominal_rejection", left, top, 100, 10000, 0, 1, COLORS[method], "7,4" if model_index else "")

    # Calibrated power at n=1000 for each model.
    for panel_index, model in enumerate(("+I+G4", "+I+R4"), start=1):
        left, top = panels[panel_index]
        axes(left, top, f"{chr(65+panel_index)}  {model} calibrated power", "n = 1,000; formula-specific training-null threshold", 0.5, 10000, 0, 1, "target branch length (log scale)", "power", False)
        for method in METHODS:
            rows = [row for row in summary if row["model"] == model and row["method"] == method and row["nsites"] == 1000]
            polyline(rows, "branch_length", "calibrated_power", left, top, 0.5, 10000, 0, 1, COLORS[method])

    # Paired calibrated power difference across n at a transition branch.
    left, top = panels[3]
    axes(left, top, "D  Paired calibrated-power difference", "refined minus legacy; t = 20", 100, 10000, -0.08, 0.08, "alignment length, n (log scale)", "power difference", True)
    zero_y = top + panel_h / 2
    svg.append(f'<line x1="{left}" x2="{left+panel_w}" y1="{zero_y:.2f}" y2="{zero_y:.2f}" stroke="#3c4043" stroke-dasharray="4,4"/>')
    for index, model in enumerate(("+I+G4", "+I+R4")):
        rows = [row for row in paired if row["model"] == model and row["rule"] == "calibrated" and row["sample"] == "alternative" and row["branch_length"] == 20]
        polyline(rows, "nsites", "refined_minus_legacy", left, top, 100, 10000, -0.08, 0.08, "#0A3FBA" if index == 0 else "#B6497D", "7,4" if index else "")

    legend = [
        ("Legacy counterfactual", COLORS["legacy_counterfactual"], ""),
        ("Refined invariant null", COLORS["refined_invariant_null"], ""),
        ("+I+G4 (solid in A; blue in D)", "#0A3FBA", ""),
        ("+I+R4 (dashed in A; pink in D)", "#B6497D", "7,4"),
    ]
    for index, (label, color, dash) in enumerate(legend):
        x = 120 + index * 275
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        svg.append(f'<line x1="{x}" x2="{x+34}" y1="88" y2="88" stroke="{color}" stroke-width="2.4"{dash_attr}/>')
        svg.append(f'<text x="{x+43}" y="92" font-size="10.5" fill="#3c4043">{label}</text>')
    svg.append("</svg>")
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    svg_path = output_prefix.with_suffix(".svg")
    svg_path.write_text("\n".join(svg) + "\n", encoding="utf-8")
    rsvg = shutil.which("rsvg-convert")
    if rsvg:
        subprocess.run([rsvg, "-f", "png", "-z", str(png_scale), "-o", str(output_prefix.with_suffix(".png")), str(svg_path)], check=True)
        subprocess.run([rsvg, "-f", "pdf", "-o", str(output_prefix.with_suffix(".pdf")), str(svg_path)], check=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--site-lengths", default="100,500,1000,5000,10000")
    parser.add_argument("--branch-lengths", default=DEFAULT_BRANCH_LENGTHS)
    parser.add_argument("--null-reps", type=int, default=6000)
    parser.add_argument("--alternative-reps", type=int, default=6000)
    parser.add_argument("--seed", type=int, default=20260711)
    parser.add_argument(
        "--outdir",
        default=str(PROJECT_ROOT / "artifacts/local/weighted-invariant-grid"),
    )
    parser.add_argument(
        "--figure-prefix",
        default=str(PROJECT_ROOT / "artifacts/local/weighted-invariant-grid/figures/overview"),
    )
    parser.add_argument("--png-scale", type=float, default=3.0)
    args = parser.parse_args()
    if args.null_reps < 100 or args.null_reps % 2:
        raise SystemExit("--null-reps must be an even integer of at least 100")

    site_lengths = [int(value) for value in args.site_lengths.split(",") if value.strip()]
    branch_lengths = [float(value) for value in args.branch_lengths.split(",") if value.strip()]
    q, pi = build_q(SUBSTITUTION_MODEL)
    eigenvalues, eigenvectors, inverse = reversible_eigendecomposition(q, pi)
    zero = int(np.argmin(np.abs(eigenvalues)))
    modes = [index for index in range(len(eigenvalues)) if index != zero]
    models = model_grid()
    summary_rows, paired_rows, assumption_rows = [], [], []
    train_count = args.null_reps // 2

    for model_name, model in models.items():
        expected_rate = float(model["proportions"] @ model["rates"])
        if abs(model["proportions"].sum() - 1.0) > 1e-12 or abs(expected_rate - 1.0) > 2e-9:
            raise AssertionError(f"Invalid {model_name} rate mixture")
        for branch_length in branch_lengths:
            null_probability, alternative_probability, scores, shift, identity_error = exact_cell(
                eigenvalues,
                eigenvectors,
                inverse,
                pi,
                modes,
                model,
                branch_length,
            )
            exact = {}
            for method in METHODS:
                mean0 = float(null_probability @ scores[method])
                variance0 = float(null_probability @ (scores[method] ** 2) - mean0 * mean0)
                mean1 = float(alternative_probability @ scores[method])
                exact[method] = (mean0, math.sqrt(variance0), mean1)
            if abs(exact["refined_invariant_null"][0]) > 2e-11:
                raise AssertionError("Refined score is not centered under the exact +I null")
            if identity_error > 3e-10:
                raise AssertionError("Refined finite/null likelihood-ratio identity failed")
            assumption_rows.append(
                {
                    "model": model_name,
                    "branch_length": branch_length,
                    "rate_mass": float(model["proportions"].sum()),
                    "mean_rate": expected_rate,
                    "probability_mass_null": float(null_probability.sum()),
                    "probability_mass_alternative": float(alternative_probability.sum()),
                    "refined_shift": shift,
                    "likelihood_identity_max_error": identity_error,
                    "refined_exact_null_mean": exact["refined_invariant_null"][0],
                    "legacy_exact_null_mean": exact["legacy_counterfactual"][0],
                }
            )

            for nsites in site_lengths:
                null_rng = np.random.default_rng(
                    seed_for(args.seed, model_name, branch_length, nsites, "null")
                )
                alternative_rng = np.random.default_rng(
                    seed_for(args.seed, model_name, branch_length, nsites, "alternative")
                )
                bootstrap_rng = np.random.default_rng(
                    seed_for(args.seed, model_name, branch_length, nsites, "bootstrap")
                )
                null_counts = simulate_counts(
                    null_probability, nsites, args.null_reps, null_rng
                )
                alternative_counts = simulate_counts(
                    alternative_probability, nsites, args.alternative_reps, alternative_rng
                )
                z_null, z_alternative, thresholds = {}, {}, {}
                for method in METHODS:
                    z_null[method] = studentized(null_counts, scores[method], nsites)
                    z_alternative[method] = studentized(
                        alternative_counts, scores[method], nsites
                    )
                    thresholds[method] = float(
                        np.quantile(z_null[method][:train_count], 0.95)
                    )
                    heldout = z_null[method][train_count:]
                    nominal_null = heldout > NOMINAL_THRESHOLD
                    calibrated_null = heldout > thresholds[method]
                    nominal_alt = z_alternative[method] > NOMINAL_THRESHOLD
                    calibrated_alt = z_alternative[method] > thresholds[method]
                    nominal_null_ci = wilson_interval(int(nominal_null.sum()), len(nominal_null))
                    calibrated_null_ci = wilson_interval(
                        int(calibrated_null.sum()), len(calibrated_null)
                    )
                    nominal_power_ci = wilson_interval(int(nominal_alt.sum()), len(nominal_alt))
                    calibrated_power_ci = wilson_interval(
                        int(calibrated_alt.sum()), len(calibrated_alt)
                    )
                    exact_mean0, exact_sd0, exact_mean1 = exact[method]
                    summary_rows.append(
                        {
                            "model": model_name,
                            "branch_length": branch_length,
                            "nsites": nsites,
                            "method": method,
                            "null_training_reps": train_count,
                            "null_heldout_reps": args.null_reps - train_count,
                            "alternative_reps": args.alternative_reps,
                            "exact_null_mean": exact_mean0,
                            "exact_null_sd": exact_sd0,
                            "exact_alternative_mean": exact_mean1,
                            "exact_standardized_separation": (exact_mean1 - exact_mean0)
                            / exact_sd0,
                            "empirical_threshold": thresholds[method],
                            "heldout_null_mean_z": float(np.mean(heldout)),
                            "heldout_null_median_z": float(np.median(heldout)),
                            "heldout_null_sd_z": float(np.std(heldout, ddof=1)),
                            "null_nominal_rejection": float(np.mean(nominal_null)),
                            "null_nominal_ci_low": nominal_null_ci[0],
                            "null_nominal_ci_high": nominal_null_ci[1],
                            "null_calibrated_rejection": float(np.mean(calibrated_null)),
                            "null_calibrated_ci_low": calibrated_null_ci[0],
                            "null_calibrated_ci_high": calibrated_null_ci[1],
                            "nominal_power": float(np.mean(nominal_alt)),
                            "nominal_power_ci_low": nominal_power_ci[0],
                            "nominal_power_ci_high": nominal_power_ci[1],
                            "calibrated_power": float(np.mean(calibrated_alt)),
                            "calibrated_power_ci_low": calibrated_power_ci[0],
                            "calibrated_power_ci_high": calibrated_power_ci[1],
                            "alternative_mean_z": float(np.mean(z_alternative[method])),
                            "alternative_median_z": float(np.median(z_alternative[method])),
                            "alternative_sd_z": float(np.std(z_alternative[method], ddof=1)),
                        }
                    )

                for rule, threshold_key in (
                    ("nominal", {method: NOMINAL_THRESHOLD for method in METHODS}),
                    ("calibrated", thresholds),
                ):
                    for sample, values in (
                        ("null_heldout", {method: z_null[method][train_count:] for method in METHODS}),
                        ("alternative", z_alternative),
                    ):
                        refined = values["refined_invariant_null"] > threshold_key["refined_invariant_null"]
                        legacy = values["legacy_counterfactual"] > threshold_key["legacy_counterfactual"]
                        low, high = paired_bootstrap_interval(
                            refined, legacy, bootstrap_rng
                        )
                        delta_z = (
                            values["refined_invariant_null"]
                            - values["legacy_counterfactual"]
                        )
                        delta_mean, delta_mean_low, delta_mean_high = mean_interval(delta_z)
                        delta_margin = (
                            values["refined_invariant_null"]
                            - threshold_key["refined_invariant_null"]
                            - values["legacy_counterfactual"]
                            + threshold_key["legacy_counterfactual"]
                        )
                        margin_mean, margin_low, margin_high = mean_interval(delta_margin)
                        paired_rows.append(
                            {
                                "model": model_name,
                                "branch_length": branch_length,
                                "nsites": nsites,
                                "sample": sample,
                                "rule": rule,
                                "replicates": len(refined),
                                "refined_fraction": float(np.mean(refined)),
                                "legacy_fraction": float(np.mean(legacy)),
                                "refined_minus_legacy": float(np.mean(refined) - np.mean(legacy)),
                                "paired_bootstrap_ci_low": low,
                                "paired_bootstrap_ci_high": high,
                                "refined_only": int(np.sum(refined & ~legacy)),
                                "legacy_only": int(np.sum(legacy & ~refined)),
                                "agreement": float(np.mean(refined == legacy)),
                                "mean_delta_z": delta_mean,
                                "mean_delta_z_ci_low": delta_mean_low,
                                "mean_delta_z_ci_high": delta_mean_high,
                                "median_delta_z": float(np.median(delta_z)),
                                "delta_z_q10": float(np.quantile(delta_z, 0.10)),
                                "delta_z_q90": float(np.quantile(delta_z, 0.90)),
                                "mean_delta_decision_margin": margin_mean,
                                "mean_delta_decision_margin_ci_low": margin_low,
                                "mean_delta_decision_margin_ci_high": margin_high,
                                "median_delta_decision_margin": float(
                                    np.median(delta_margin)
                                ),
                            }
                        )

    outdir = Path(args.outdir)
    write_tsv(outdir / "weighted_invariant_summary.tsv", summary_rows)
    write_tsv(outdir / "weighted_invariant_paired.tsv", paired_rows)
    write_tsv(outdir / "weighted_invariant_assumptions.tsv", assumption_rows)
    render_figure(
        summary_rows,
        paired_rows,
        Path(args.figure_prefix),
        args.png_scale,
    )
    print(outdir / "weighted_invariant_summary.tsv")
    print(outdir / "weighted_invariant_paired.tsv")
    print(outdir / "weighted_invariant_assumptions.tsv")


if __name__ == "__main__":
    main()
