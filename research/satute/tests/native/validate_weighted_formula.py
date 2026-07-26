#!/usr/bin/env python3
"""Validate the eigenvalue_weighted SatuTe formula against dominant.

Consolidates the experiment chain that diagnosed and confirmed the fix for
tree/satute/satute_statistics.cpp's centring change: the raw cross-moment
statistic sum_k w_k L_k R_k is only a valid independence test when
E[L_k]=E[R_k]=0, which fails under a misspecified substitution model and
biases the statistic in a way that GROWS with alignment length. The fix
centres the statistic into a covariance, which is zero under independence
regardless of the marginal means.

Every experiment below is independently runnable (--which) and self-checking
(asserts its own expected outcome; exits non-zero on regression). Read
EXPERIMENTS at the bottom for the full list and what each one establishes.

Usage:
    research/satute/satute developer validate-weighted-formula -- \\
        --iqtree ./build/iqtree3 --outdir /tmp/satute-weighted-validate

    # single experiment, fast iteration:
    ... --which decoupled_weight_calibration --reps-scale 0.2

Design notes shared across the simulation-based experiments:
  * The tested branch's own length never enters the likelihood computation
    (computeOrdinaryPartial recurses away from it via left_dad/right_dad) --
    only buildFormulaWeights reads it. This is what makes the "decoupled"
    design possible: simulate the tested branch at length T_SIM (independent
    sides, true null) then re-analyze the SAME alignment via `-te -blfix`
    declaring a DIFFERENT length T_weight, so the formula's weight vector is
    whatever we choose while the true correlation stays exactly what the
    simulation produced.
  * TRUE_MODEL = GTR{1,2,1,1,2,1}+F{0.30,0.20,0.20,0.30} is mathematically
    HKY85 (transitions rate 2, transversions rate 1, skewed frequencies) with
    a repeated secondary eigenvalue (2 nonzero eigenvalues: one dominant,
    one degenerate pair). K2P is the natural misspecification of this model
    -- same rate-ratio structure, but forces equal (0.25 each) frequencies.
  * Experiments 1-8 all run at that ONE point in eigen-structure space
    (gap -0.51, singleton dominant mode). Experiment 9 is the only one that
    varies the spectrum itself; see its docstring for why that matters.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import math
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))
sys.path.insert(0, str(PROJECT_ROOT / "tests" / "exact"))

from satute_analysis.models import build_q, reversible_eigendecomposition  # noqa: E402
from test_saturation_null import exact_patterns  # noqa: E402


TRUE_MODEL = "GTR{1,2,1,1,2,1}+F{0.30,0.20,0.20,0.30}"
MISSPEC_MODEL = "K2P"
ALPHA = 0.05
SPLIT = "A,B"


# --------------------------------------------------------------------------
# Shared infrastructure
# --------------------------------------------------------------------------

def stable_seed(*parts: object) -> int:
    """Deterministic seed from arbitrary parts, so experiments composed in
    any order or run individually never collide (unlike hand-picked additive
    offsets, which is what the ad hoc predecessor scripts used)."""
    digest = hashlib.sha256("|".join(str(p) for p in parts).encode()).digest()
    return int.from_bytes(digest[:4], "big") % 900_000_000 + 1


def tree_newick(a: float, b: float, c: float, d: float, e: float, t: float) -> str:
    """The 5-taxon design used throughout: split A,B | C,D,E across the
    tested branch of length t, with external branches a,b,c,d,e."""
    return f"(A:{a:.6g},B:{b:.6g},(C:{c:.6g},D:{d:.6g},E:{e:.6g}):{t:.6g});\n"


def run(cmd: list[str]) -> None:
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"command failed ({result.returncode}): {' '.join(cmd)}\n{result.stdout}")


def simulate_alignment(iqtree: str, workdir: Path, tree_path: Path, model: str, nsites: int, seed: int) -> Path:
    prefix = workdir / "sim"
    run([iqtree, "--alisim", str(prefix), "-t", str(tree_path), "-m", model, "--length", str(nsites),
         "--seed", str(seed), "-af", "fasta", "--quiet", "--skip-bl-check"])
    return Path(str(prefix) + ".fa")


def run_satute(iqtree: str, alignment: Path, tree_path: Path, model: str, seed: int, prefix: Path) -> Path:
    run([iqtree, "-s", str(alignment), "-te", str(tree_path), "-m", model, "-blfix", "--satute",
         "--prefix", str(prefix), "-T", "1", "--redo", "--quiet", "--seed", str(seed)])
    return Path(str(prefix) + ".sat.stat")


def parse_pooled_pvalues(stat_path: Path, split: str = SPLIT) -> dict[str, float]:
    """{formula: satP} for the pooled row on the requested split. Raises if
    the file is missing or unparsable rather than silently returning {} --
    the ad hoc predecessor scripts swallowed this, which can hide a run that
    silently produced no rows at all."""
    header = None
    out: dict[str, float] = {}
    for line in stat_path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if fields[0] == "ID":
            header = fields
            continue
        if header is None:
            raise ValueError(f"{stat_path}: data row before header")
        row = dict(zip(header, fields))
        if row.get("Split") == split and row.get("RateCategory") == "pooled":
            out[row["Formula"]] = float(row["satP"])
    return out


def parse_pooled_rows(stat_path: Path, split: str = SPLIT) -> dict[str, dict[str, str]]:
    """{formula: whole row} for the pooled row on the requested split.

    Needed where the p-value alone is not enough -- notably the ``Modes``
    column, which records WHICH eigen-modes the formula actually selected.
    """
    header = None
    out: dict[str, dict[str, str]] = {}
    for line in stat_path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if fields[0] == "ID":
            header = fields
            continue
        if header is None:
            raise ValueError(f"{stat_path}: data row before header")
        row = dict(zip(header, fields))
        if row.get("Split") == split and row.get("RateCategory") == "pooled":
            out[row["Formula"]] = row
    return out


def rejection_ci(p_values: np.ndarray, alpha: float = ALPHA) -> tuple[float, float, float, int]:
    """(rejection_rate, ci_low, ci_high, n) using a normal-approximation 95% CI
    around the NOMINAL rate -- appropriate here because we are asking "is the
    observed rate distinguishable from nominal", not estimating an unknown rate."""
    n = len(p_values)
    rate = float((p_values <= alpha).mean())
    se = 1.96 * math.sqrt(alpha * (1 - alpha) / n)
    # Clip for display/sanity: a rejection rate cannot be negative or exceed
    # 1, but the normal approximation can cross those bounds at small n.
    return rate, max(0.0, alpha - se), min(1.0, alpha + se), n


def mcnemar(a_rejects: np.ndarray, b_rejects: np.ndarray) -> tuple[int, int, float]:
    """Exact paired McNemar test via the binomial: b_only = a rejects, b does
    not; c_only = b rejects, a does not. Returns (b_only, c_only, two-sided p)."""
    b_only = int((a_rejects & ~b_rejects).sum())
    c_only = int((~a_rejects & b_rejects).sum())
    total = b_only + c_only
    pvalue = stats.binomtest(c_only, total, 0.5).pvalue if total > 0 else 1.0
    return b_only, c_only, pvalue


def parallel_map(fn, items, workers: int | None = None) -> list:
    with ThreadPoolExecutor(max_workers=workers or os.cpu_count()) as executor:
        return list(executor.map(fn, items))


def eigen_gap(model: str = TRUE_MODEL) -> tuple[float, float]:
    """(dominant eigenvalue, gap to the next-largest nonzero eigenvalue)."""
    q, pi = build_q(model)
    eigenvalues, _, _ = reversible_eigendecomposition(q, pi)
    nonzero = np.sort(eigenvalues[np.abs(eigenvalues) > 1e-9])[::-1]
    return float(nonzero[0]), float(nonzero[1] - nonzero[0])


def section(title: str) -> None:
    print(f"\n{'=' * len(title)}\n{title}\n{'=' * len(title)}")


# --------------------------------------------------------------------------
# Experiment 1 -- exact per-site score skewness under the saturated null
# --------------------------------------------------------------------------

def exact_skewness(args: argparse.Namespace) -> dict:
    """Was the small GTR anti-conservatism seen pre-fix explained by skewness
    of the per-site score (slow normal convergence)? Computed exactly by
    pattern enumeration, not simulated. Expected: skewness ~ 0 (this
    hypothesis was ruled out; see gtr_null_scaling for what it actually was)."""
    section("1. Exact per-site score skewness under H0")
    p, _cat, left, right, _soft = exact_patterns()
    x = left * right
    mean = float(p @ x)
    var = float(p @ (x - mean) ** 2)
    skew = float(p @ (x - mean) ** 3) / var**1.5
    kurt = float(p @ (x - mean) ** 4) / var**2 - 3.0
    print(f"  mean = {mean:.3e} (H0 requires 0)")
    print(f"  skewness = {skew:+.4f}, excess kurtosis = {kurt:+.3f}")
    ok = abs(skew) < 1e-6 and abs(mean) < 1e-9
    print(f"  {'PASS' if ok else 'FAIL'}: skewness negligible -> not the source of any observed miscalibration")
    assert ok, f"unexpected nonzero skew/mean: skew={skew}, mean={mean}"
    return {"mean": mean, "skewness": skew, "excess_kurtosis": kurt}


# --------------------------------------------------------------------------
# Experiment 2 -- exact variance-estimator bias
# --------------------------------------------------------------------------

def exact_variance_bias(args: argparse.Namespace) -> dict:
    """SatuTe's variance estimator factorizes Var(L*R) as E[L^2]*E[R^2],
    assuming L independent of R. Exact size of that approximation error,
    computed by pattern enumeration. Expected: a fraction of a percent --
    negligible, ruling out "replace the variance formula" as a fix for the
    misspecification bias (which is a bias in the MEAN, not the variance;
    see misspecification_null)."""
    section("2. Exact variance-estimator bias (factorized vs true)")
    p, _cat, left, right, _soft = exact_patterns()
    x = left * right
    true_var = float(p @ (x**2)) - float(p @ x) ** 2
    plug_in = float(p @ (left**2)) * float(p @ (right**2))
    ratio = plug_in / true_var
    print(f"  true Var(L*R)              = {true_var:.6e}")
    print(f"  estimator target E[L^2]E[R^2] = {plug_in:.6e}")
    print(f"  ratio = {ratio:.4f} ({100 * abs(1 - ratio):.2f}% {'under' if ratio < 1 else 'over'}-statement)")
    ok = abs(1 - ratio) < 0.02
    print(f"  {'PASS' if ok else 'FAIL'}: bias under 2%, i.e. negligible relative to sampling noise at any n")
    assert ok, f"variance-estimator bias larger than expected: ratio={ratio}"
    return {"true_variance": true_var, "estimator_target": plug_in, "ratio": ratio}


# --------------------------------------------------------------------------
# Experiment 3 -- correctly-specified GTR null, does deviation shrink with n
# --------------------------------------------------------------------------

def gtr_null_scaling(args: argparse.Namespace) -> dict:
    """Under the CORRECTLY specified model, is residual anti-conservatism at
    finite n a shrinking finite-sample artifact (fine) or does it persist
    (would indicate a real defect)? This is the control for
    misspecification_scaling below, which asks the same question under a
    WRONG model and gets the opposite answer."""
    section("3. Correctly-specified GTR null: does deviation shrink with n?")
    reps = max(30, int(200 * args.reps_scale))
    site_counts = (250, 1000, 4000, 16000)
    workdir = args.outdir / "gtr_null_scaling"
    workdir.mkdir(parents=True, exist_ok=True)

    def one(task):
        nsites, rep = task
        d = workdir / f"n{nsites}" / f"r{rep}"
        d.mkdir(parents=True, exist_ok=True)
        tree = d / "t.nwk"
        tree.write_text(tree_newick(0.05, 0.05, 0.05, 0.05, 0.05, 100.0))
        seed = stable_seed("gtr_null_scaling", nsites, rep)
        alignment = simulate_alignment(args.iqtree, d, tree, TRUE_MODEL, nsites, seed)
        stat = run_satute(args.iqtree, alignment, tree, TRUE_MODEL, seed, d / "s")
        return nsites, parse_pooled_pvalues(stat)

    tasks = [(n, r) for n in site_counts for r in range(1, reps + 1)]
    results: dict[int, list[float]] = {n: [] for n in site_counts}
    for nsites, pvals in parallel_map(one, tasks, args.workers):
        if "eigenvalue_weighted" in pvals:
            results[nsites].append(pvals["eigenvalue_weighted"])

    print(f"{'sites':>7} {'reject':>8} {'95% CI':>18}")
    rates = {}
    for nsites in site_counts:
        rate, lo, hi, n = rejection_ci(np.array(results[nsites]))
        rates[nsites] = rate
        flag = "" if lo <= rate <= hi else "  <-- outside"
        print(f"{nsites:>7} {rate:>8.4f} {f'[{lo:.4f},{hi:.4f}]':>18}{flag}   n={n}")

    # The deviation at the smallest n should not be larger at the largest n.
    ok = abs(rates[site_counts[-1]] - ALPHA) <= abs(rates[site_counts[0]] - ALPHA) + 0.02
    print(f"  {'PASS' if ok else 'FAIL'}: deviation does not grow with n (finite-sample, correctly specified)")
    assert ok, f"deviation grew with n under the CORRECT model: {rates}"
    return {"rejection_by_nsites": rates}


# --------------------------------------------------------------------------
# Experiment 4 -- calibration across the weight spectrum (correct model)
# --------------------------------------------------------------------------

def decoupled_weight_calibration(args: argparse.Namespace) -> dict:
    """Is eigenvalue_weighted correctly calibrated when its own weight vector
    is genuinely non-degenerate, under a correctly specified model? Uses the
    decoupling design (see module docstring): true independence from
    T_SIM=100, weight varied via T_weight in the analysis tree. dominant
    ignores branch length entirely (weight is always 1), so its rejection
    rate must be IDENTICAL across every T_weight -- that is the built-in
    sanity check that only the weight actually changed."""
    section("4. Calibration across the weight spectrum (correctly specified)")
    reps = max(200, int(2000 * args.reps_scale))
    t_sim = 100.0
    t_weights = (0.3, 1.0, 3.0, 100.0)
    workdir = args.outdir / "decoupled_weight_calibration"
    workdir.mkdir(parents=True, exist_ok=True)

    def simulate_rep(rep):
        d = workdir / f"r{rep}"
        d.mkdir(parents=True, exist_ok=True)
        sim_tree = d / "sim.nwk"
        sim_tree.write_text(tree_newick(0.05, 0.05, 0.05, 0.05, 0.05, t_sim))
        seed = stable_seed("decoupled_weight_calibration", rep)
        alignment = simulate_alignment(args.iqtree, d, sim_tree, TRUE_MODEL, 1000, seed)
        return d, alignment, seed

    def analyze(args_tuple):
        d, alignment, seed, t_w = args_tuple
        analysis_tree = d / f"analysis_T{t_w:g}.nwk"
        analysis_tree.write_text(tree_newick(0.05, 0.05, 0.05, 0.05, 0.05, t_w))
        stat = run_satute(args.iqtree, alignment, analysis_tree, TRUE_MODEL, seed, d / f"s_T{t_w:g}")
        return t_w, parse_pooled_pvalues(stat)

    sims = parallel_map(simulate_rep, range(1, reps + 1), args.workers)
    tasks = [(d, alignment, seed, t_w) for (d, alignment, seed) in sims for t_w in t_weights]
    acc: dict[float, dict[str, list[float]]] = {t_w: {"dominant": [], "eigenvalue_weighted": []} for t_w in t_weights}
    for t_w, pvals in parallel_map(analyze, tasks, args.workers):
        for formula in ("dominant", "eigenvalue_weighted"):
            if formula in pvals:
                acc[t_w][formula].append(pvals[formula])

    dominant_lambda, gap = eigen_gap()
    print(f"  dominant eigenvalue={dominant_lambda:.4f}  gap to secondary={gap:+.4f}\n")
    print(f"{'T_weight':>9} {'weight(secondary)':>18} {'dominant reject':>16} {'weighted reject':>16} {'95% CI':>16}")
    dominant_rates = set()
    all_ok = True
    for t_w in t_weights:
        weight = math.exp(gap * t_w)
        dominant = np.array(acc[t_w]["dominant"])
        weighted = np.array(acc[t_w]["eigenvalue_weighted"])
        rd = float((dominant <= ALPHA).mean())
        rw, lo, hi, _n = rejection_ci(weighted)
        dominant_rates.add(round(rd, 4))
        ok = lo <= rw <= hi
        all_ok &= ok
        flag = "" if ok else "  <-- outside"
        print(f"{t_w:>9g} {weight:>18.4f} {rd:>16.4f} {rw:>16.4f} {f'[{lo:.4f},{hi:.4f}]':>16}{flag}")

    dominant_invariant = len(dominant_rates) == 1
    print(f"\n  sanity check -- dominant identical across all T_weight: {'PASS' if dominant_invariant else 'FAIL'} {sorted(dominant_rates)}")
    print(f"  {'PASS' if all_ok else 'FAIL'}: weighted calibrated at every weight level")
    assert dominant_invariant, f"dominant should be invariant to T_weight, got {dominant_rates}"
    assert all_ok, "weighted formula miscalibrated at some non-degenerate weight"
    return {"dominant_rate": next(iter(dominant_rates)), "weighted_by_T_weight": {t: float((np.array(acc[t]['eigenvalue_weighted']) <= ALPHA).mean()) for t in t_weights}}


# --------------------------------------------------------------------------
# Experiment 5 -- matched-branch power, with significance testing
# --------------------------------------------------------------------------

def matched_power(args: argparse.Namespace) -> dict:
    """Does eigenvalue_weighted have more power than dominant when there is
    real signal, at matched (non-decoupled) branch lengths spanning the
    informative-to-saturated transition? Paired per-replicate (same
    alignment analyzed by both formulas), so McNemar's exact test applies --
    this is the check the ad hoc predecessor script (power_gtr.py) never
    ran, and the nominal edge it reported did not survive it (n=500:
    p=0.21-0.26, not significant)."""
    section("5. Matched-branch power: dominant vs weighted, with significance")
    reps = max(100, int(500 * args.reps_scale))
    branches = (0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0)
    workdir = args.outdir / "matched_power"
    workdir.mkdir(parents=True, exist_ok=True)

    def one(task):
        t, rep = task
        d = workdir / f"T{t:g}" / f"r{rep}"
        d.mkdir(parents=True, exist_ok=True)
        tree = d / "t.nwk"
        tree.write_text(tree_newick(0.05, 0.05, 0.05, 0.05, 0.05, t))
        seed = stable_seed("matched_power", t, rep)
        alignment = simulate_alignment(args.iqtree, d, tree, TRUE_MODEL, 1000, seed)
        stat = run_satute(args.iqtree, alignment, tree, TRUE_MODEL, seed, d / "s")
        return t, parse_pooled_pvalues(stat)

    tasks = [(t, r) for t in branches for r in range(1, reps + 1)]
    acc: dict[float, dict[str, list[float]]] = {t: {"dominant": [], "eigenvalue_weighted": []} for t in branches}
    for t, pvals in parallel_map(one, tasks, args.workers):
        for formula in ("dominant", "eigenvalue_weighted"):
            if formula in pvals:
                acc[t][formula].append(pvals[formula])

    _dominant_lambda, gap = eigen_gap()
    print(f"{'branch':>7} {'weight':>7} {'dominant':>9} {'weighted':>9} {'delta':>7} {'McNemar p':>10}  verdict")
    report = {}
    for t in branches:
        dominant = (np.array(acc[t]["dominant"]) <= ALPHA)
        weighted = (np.array(acc[t]["eigenvalue_weighted"]) <= ALPHA)
        rd, rw = float(dominant.mean()), float(weighted.mean())
        b_only, c_only, pvalue = mcnemar(dominant, weighted)
        verdict = "significant" if pvalue < ALPHA else "not significant"
        report[t] = {"dominant": rd, "weighted": rw, "delta": rw - rd, "mcnemar_p": pvalue,
                     "dominant_only": b_only, "weighted_only": c_only}
        print(f"{t:>7g} {math.exp(gap * t):>7.3f} {rd:>9.3f} {rw:>9.3f} {rw - rd:>+7.3f} {pvalue:>10.4f}  {verdict}")
    print("\n  This experiment reports, it does not assert -- whether weighted has\n"
          "  a real power edge is an open empirical question, not a fixed fact.")
    return report


# --------------------------------------------------------------------------
# Experiment 6 -- misspecification: the decisive decoupled test
# --------------------------------------------------------------------------

def misspecification_null(args: argparse.Namespace) -> dict:
    """The test that actually distinguishes "more power" from "inflated Type
    I error" under misspecification: true correlation forced to exactly zero
    (decoupled design) + non-degenerate weight + K2P instead of the true
    model. If eigenvalue_weighted's rejection rate departs from nominal here,
    it is pure bias -- there is no signal in this data by construction.
    Paired per replicate (same alignment, both models) for a clean z-shift
    estimate."""
    section("6. Misspecification, decoupled null (true independence, K2P vs correct)")
    reps = max(200, int(1500 * args.reps_scale))
    t_sim, t_weight = 100.0, 3.0
    workdir = args.outdir / "misspecification_null"
    workdir.mkdir(parents=True, exist_ok=True)

    def one(rep):
        d = workdir / f"r{rep}"
        d.mkdir(parents=True, exist_ok=True)
        sim_tree = d / "sim.nwk"
        sim_tree.write_text(tree_newick(0.05, 0.05, 0.05, 0.05, 0.05, t_sim))
        seed = stable_seed("misspecification_null", rep)
        alignment = simulate_alignment(args.iqtree, d, sim_tree, TRUE_MODEL, 1000, seed)
        analysis_tree = d / "analysis.nwk"
        analysis_tree.write_text(tree_newick(0.05, 0.05, 0.05, 0.05, 0.05, t_weight))
        results = {}
        for arm, model in (("correct", TRUE_MODEL), ("misspec", MISSPEC_MODEL)):
            stat = run_satute(args.iqtree, alignment, analysis_tree, model, seed, d / f"s_{arm}")
            results[arm] = parse_pooled_pvalues(stat)
        return results

    all_results = parallel_map(one, range(1, reps + 1), args.workers)

    report = {}
    print(f"{'arm':>10} {'formula':>20} {'reject':>8} {'95% CI':>18}")
    for arm in ("correct", "misspec"):
        for formula in ("dominant", "eigenvalue_weighted"):
            p = np.array([r[arm][formula] for r in all_results if formula in r[arm]])
            rate, lo, hi, n = rejection_ci(p)
            report[f"{arm}_{formula}"] = rate
            flag = "" if lo <= rate <= hi else "  <-- OUTSIDE"
            print(f"{arm:>10} {formula:>20} {rate:>8.4f} {f'[{lo:.4f},{hi:.4f}]':>18}{flag}")

    correct_p = np.array([r["correct"]["eigenvalue_weighted"] for r in all_results])
    misspec_p = np.array([r["misspec"]["eigenvalue_weighted"] for r in all_results])
    z_correct = stats.norm.isf(correct_p)
    z_misspec = stats.norm.isf(misspec_p)
    shift = float(z_misspec.mean() - z_correct.mean())
    print(f"\n  mean Z (weighted): correct={z_correct.mean():+.4f}  misspec={z_misspec.mean():+.4f}  shift={shift:+.4f}")

    _, lo, hi, _n = rejection_ci(misspec_p)
    ok = lo <= report["misspec_eigenvalue_weighted"] <= hi
    print(f"  {'PASS' if ok else 'FAIL'}: misspecified weighted stays inside the null CI (post-centring-fix)")
    assert ok, f"misspecification bias regressed: rate={report['misspec_eigenvalue_weighted']}, CI=[{lo},{hi}]"
    report["z_shift"] = shift
    return report


# --------------------------------------------------------------------------
# Experiment 7 -- does the misspecification bias grow or shrink with n
# --------------------------------------------------------------------------

def misspecification_scaling(args: argparse.Namespace) -> dict:
    """Companion to gtr_null_scaling, under the WRONG model. Pre-fix this
    grew with n (250->16000 sites: 0.064->0.091->0.178->0.380 -- structural,
    not finite-sample). Post-fix it should stay flat, matching the
    correctly-specified control."""
    section("7. Misspecification bias vs alignment length (structural, not finite-sample)")
    reps = max(200, int(800 * args.reps_scale))
    t_sim, t_weight = 100.0, 3.0
    site_counts = (250, 1000, 4000, 16000)
    workdir = args.outdir / "misspecification_scaling"
    workdir.mkdir(parents=True, exist_ok=True)

    def one(task):
        nsites, rep = task
        d = workdir / f"n{nsites}" / f"r{rep}"
        d.mkdir(parents=True, exist_ok=True)
        sim_tree = d / "sim.nwk"
        sim_tree.write_text(tree_newick(0.05, 0.05, 0.05, 0.05, 0.05, t_sim))
        seed = stable_seed("misspecification_scaling", nsites, rep)
        alignment = simulate_alignment(args.iqtree, d, sim_tree, TRUE_MODEL, nsites, seed)
        analysis_tree = d / "analysis.nwk"
        analysis_tree.write_text(tree_newick(0.05, 0.05, 0.05, 0.05, 0.05, t_weight))
        stat = run_satute(args.iqtree, alignment, analysis_tree, MISSPEC_MODEL, seed, d / "s")
        return nsites, parse_pooled_pvalues(stat)

    tasks = [(n, r) for n in site_counts for r in range(1, reps + 1)]
    results: dict[int, list[float]] = {n: [] for n in site_counts}
    for nsites, pvals in parallel_map(one, tasks, args.workers):
        if "eigenvalue_weighted" in pvals:
            results[nsites].append(pvals["eigenvalue_weighted"])

    print(f"{'sites':>7} {'reject':>8} {'95% CI':>18}")
    rates = {}
    for nsites in site_counts:
        rate, lo, hi, n = rejection_ci(np.array(results[nsites]))
        rates[nsites] = rate
        flag = "" if lo <= rate <= hi else "  <-- outside"
        print(f"{nsites:>7} {rate:>8.4f} {f'[{lo:.4f},{hi:.4f}]':>18}{flag}   n={n}")

    ok = abs(rates[site_counts[-1]] - ALPHA) <= abs(rates[site_counts[0]] - ALPHA) + 0.03
    print(f"  {'PASS' if ok else 'FAIL'}: deviation does not grow with n under misspecification (post-fix)")
    assert ok, f"misspecification bias grows with n -- regression: {rates}"
    return {"rejection_by_nsites": rates}


# --------------------------------------------------------------------------
# Experiment 8 -- IUPAC ambiguity-code fix, unrelated to the weighted-formula
# work above but validated in the same session and belongs in the same suite
# --------------------------------------------------------------------------

def ambiguity_fix(args: argparse.Namespace) -> dict:
    """Does expanding IUPAC ambiguity codes (computeTipLikelihood) move
    results toward ground truth, versus the old behaviour of treating them
    as fully unknown? Same alignment three ways (truth / R-masked /
    N-masked), identical fixed tree, so the only thing that varies is how
    leaf partials treat the masked sites."""
    section("8. IUPAC ambiguity-code fix: closer to ground truth?")
    reps = max(10, int(60 * args.reps_scale))
    mask_taxon, frac = "C", 0.40
    tree_str = tree_newick(0.05, 0.05, 0.15, 0.05, 0.05, 0.4)
    workdir = args.outdir / "ambiguity_fix"
    workdir.mkdir(parents=True, exist_ok=True)

    def read_fasta(p: Path) -> dict[str, str]:
        seqs: dict[str, str] = {}
        name = None
        for line in p.read_text().splitlines():
            if line.startswith(">"):
                name = line[1:].strip()
                seqs[name] = ""
            elif name:
                seqs[name] += line.strip()
        return seqs

    def write_fasta(seqs: dict[str, str], p: Path) -> None:
        p.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))

    def satz(fasta: Path, tree: Path, seed: int, prefix: Path) -> float:
        stat = run_satute(args.iqtree, fasta, tree, "JC", seed, prefix)
        header, z_value = None, None
        for line in stat.read_text().splitlines():
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[0] == "ID":
                header = fields
                continue
            row = dict(zip(header, fields))
            if row.get("Split") == SPLIT and row.get("Formula") == "dominant" and row.get("RateCategory") == "pooled":
                z_value = float(row["satZ"])
        return z_value

    def one(rep):
        import random
        d = workdir / f"r{rep}"
        d.mkdir(parents=True, exist_ok=True)
        tree = d / "t.nwk"
        tree.write_text(tree_str)
        seed = stable_seed("ambiguity_fix", rep)
        alignment = simulate_alignment(args.iqtree, d, tree, "JC", 2000, seed)
        seqs = read_fasta(alignment)
        rng = random.Random(seed)
        target = list(seqs[mask_taxon])
        purine_pos = [i for i, ch in enumerate(target) if ch in "AG"]
        chosen = rng.sample(purine_pos, int(len(purine_pos) * frac))
        for label, replacement in (("R", "R"), ("N", "N")):
            masked = dict(seqs)
            s = list(target)
            for i in chosen:
                s[i] = replacement
            masked[mask_taxon] = "".join(s)
            write_fasta(masked, d / f"{label}.fa")
        z_true = satz(alignment, tree, seed, d / "truth")
        z_r = satz(d / "R.fa", tree, seed, d / "amb_R")
        z_n = satz(d / "N.fa", tree, seed, d / "amb_N")
        return abs(z_r - z_true), abs(z_n - z_true)

    deltas = parallel_map(one, range(1, reps + 1), args.workers)
    delta_r = np.array([d[0] for d in deltas])
    delta_n = np.array([d[1] for d in deltas])
    closer = int((delta_r < delta_n).sum())
    print(f"  replicates={len(delta_r)}  masked {frac:.0%} of taxon {mask_taxon}'s purine sites")
    print(f"  mean |satZ-truth|, R (expanded to A/G): {delta_r.mean():.4f}")
    print(f"  mean |satZ-truth|, N (fully unknown):   {delta_n.mean():.4f}")
    print(f"  R closer to truth in {closer}/{len(delta_r)} replicates")
    ok = closer > len(delta_r) * 0.7 and delta_r.mean() < delta_n.mean()
    print(f"  {'PASS' if ok else 'FAIL'}")
    assert ok, "ambiguity fix no longer closer to ground truth than the old behaviour"
    return {"mean_abs_error_R": float(delta_r.mean()), "mean_abs_error_N": float(delta_n.mean()), "closer_fraction": closer / len(delta_r)}


# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# Experiment 9 -- spectral gap: coverage of the eigen-structure, and the
# mode-selection stability of dominant vs eigenvalue_weighted
# --------------------------------------------------------------------------

def spectral_gap_stability(args: argparse.Namespace) -> dict:
    """Every other experiment here runs at ONE point in eigen-structure space
    (TRUE_MODEL: gap -0.51, a singleton dominant mode). That leaves the rest
    of the space untested, including a whole branch of findDominantModes.

    Part A sweeps kappa in GTR{1,kappa,1,1,kappa,1} with equal frequencies.
    This one parameter moves the spectral gap continuously through zero AND
    changes the dominant mode's multiplicity, because the spectrum is a
    singleton plus a degenerate pair:
        kappa < 1  -> the PAIR is dominant   (multiplicity 2)
        kappa = 1  -> JC, all degenerate     (multiplicity 3)
        kappa > 1  -> the SINGLETON dominant (multiplicity 1)

    Part B is the reason this matters. findDominantModes selects modes within
    a hard tolerance max(1e-8, |lambda|*1e-6) of the largest non-zero
    eigenvalue -- a DISCONTINUOUS choice. When the analysis model's kappa is
    ML-estimated and the truth is near-degenerate, kappa_hat lands on a random
    side of 1, so `dominant` picks a different mode set per replicate and its
    estimand changes between runs of the same experiment. eigenvalue_weighted
    has no selection step (it always takes every non-zero mode, with
    continuously varying weights) and so cannot exhibit this.

    Part B asserts only that eigenvalue_weighted is stable. `dominant`'s
    instability is REPORTED, not asserted: it is a property of the current
    published formula, and a future fix that removes it must not fail this
    suite.
    """
    section("9. Spectral gap: eigen-structure coverage and mode-selection stability")

    # ---- Part A: analytic sweep, no simulation needed ---------------------
    print("Part A -- reachable eigen-structure (analytic)\n")
    print(f"{'kappa':>7} {'nonzero eigenvalues':>34} {'gap':>9} {'dominant modes':>15}")
    multiplicities = {}
    for kappa in (0.2, 0.5, 0.95, 1.0, 1.05, 2.0, 5.0):
        model = f"GTR{{1,{kappa},1,1,{kappa},1}}+F{{0.25,0.25,0.25,0.25}}"
        q, pi = build_q(model)
        eigenvalues, _, _ = reversible_eigendecomposition(q, pi)
        nonzero = np.sort(eigenvalues[np.abs(eigenvalues) > 1e-9])[::-1]
        lam = float(nonzero.max())
        tol = max(1e-8, abs(lam) * 1e-6)
        n_dominant = int(np.sum(np.abs(nonzero - lam) <= tol))
        gap = float(np.sort(nonzero)[::-1][n_dominant] - lam) if n_dominant < len(nonzero) else 0.0
        multiplicities[kappa] = n_dominant
        print(f"{kappa:>7g} {str(np.round(nonzero, 5)):>34} {gap:>9.5f} {n_dominant:>15}")

    covered = set(multiplicities.values())
    ok_coverage = {1, 2, 3} <= covered
    print(f"\n  multiplicities covered by this sweep: {sorted(covered)}")
    print(f"  {'PASS' if ok_coverage else 'FAIL'}: sweep exercises singleton, degenerate-pair, "
          "and fully-degenerate dominant sets")
    assert ok_coverage, f"kappa sweep failed to cover multiplicities 1/2/3, got {sorted(covered)}"

    # ---- Part B: stability under an ML-estimated, structurally degenerate model
    reps = max(24, int(80 * args.reps_scale))
    print(f"\nPart B -- mode-set stability, JC truth analysed under free {MISSPEC_MODEL} "
          f"(kappa estimated), {reps} replicates\n")
    workdir = args.outdir / "spectral_gap_stability"
    workdir.mkdir(parents=True, exist_ok=True)

    def one(rep):
        d = workdir / f"r{rep}"
        d.mkdir(parents=True, exist_ok=True)
        tree = d / "t.nwk"
        tree.write_text(tree_newick(0.05, 0.05, 0.05, 0.05, 0.05, 1.0))
        seed = stable_seed("spectral_gap_stability", rep)
        alignment = simulate_alignment(args.iqtree, d, tree, "JC", 1000, seed)
        stat = run_satute(args.iqtree, alignment, tree, MISSPEC_MODEL, seed, d / "s")
        return parse_pooled_rows(stat)

    rows = parallel_map(one, range(1, reps + 1), args.workers)

    report: dict = {"multiplicity_by_formula": {}}
    for formula in ("dominant", "eigenvalue_weighted"):
        by_multiplicity: dict[int, list[float]] = {}
        for row in rows:
            if formula not in row:
                continue
            n_modes = len(row[formula]["Modes"].split(","))
            by_multiplicity.setdefault(n_modes, []).append(float(row[formula]["satZ"]))
        counts = {k: len(v) for k, v in sorted(by_multiplicity.items())}
        pooled = np.array([z for values in by_multiplicity.values() for z in values])
        print(f"  {formula:>20}: multiplicity {counts}   "
              f"satZ mean={pooled.mean():+.3f} sd={pooled.std():.3f}")
        report["multiplicity_by_formula"][formula] = counts

        if len(by_multiplicity) > 1:
            within = float(np.mean([np.std(v) for v in by_multiplicity.values() if len(v) > 1]))
            inflation = float(pooled.std() / within) if within > 0 else float("nan")
            for n_modes, values in sorted(by_multiplicity.items()):
                arr = np.array(values)
                print(f"      {n_modes}-mode subset (n={len(arr):>3}): "
                      f"satZ mean={arr.mean():+.3f} sd={arr.std():.3f}")
            print(f"      -> mode-set flip inflates satZ sd by {inflation:.2f}x "
                  f"(within-group {within:.3f} vs pooled {pooled.std():.3f})")
            report[f"{formula}_variance_inflation"] = inflation

    weighted_counts = report["multiplicity_by_formula"].get("eigenvalue_weighted", {})
    weighted_stable = len(weighted_counts) == 1
    dominant_counts = report["multiplicity_by_formula"].get("dominant", {})
    print(f"\n  {'PASS' if weighted_stable else 'FAIL'}: eigenvalue_weighted selects one fixed "
          f"mode set across all replicates {weighted_counts}")
    if len(dominant_counts) > 1:
        print(f"  NOTE (reported, not asserted): dominant's mode set is unstable "
              f"{dominant_counts} -- its estimand depends on estimation noise in kappa_hat.")
    assert weighted_stable, (
        f"eigenvalue_weighted mode set is no longer stable: {weighted_counts}"
    )
    return report


EXPERIMENTS = {
    "exact_skewness": exact_skewness,
    "exact_variance_bias": exact_variance_bias,
    "gtr_null_scaling": gtr_null_scaling,
    "decoupled_weight_calibration": decoupled_weight_calibration,
    "matched_power": matched_power,
    "misspecification_null": misspecification_null,
    "misspecification_scaling": misspecification_scaling,
    "ambiguity_fix": ambiguity_fix,
    "spectral_gap_stability": spectral_gap_stability,
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--iqtree", required=True, help="Path to the iqtree3 binary.")
    parser.add_argument("--outdir", required=True, type=Path, help="Working/output directory.")
    parser.add_argument("--which", default="all",
                         help="Comma-separated experiment names, or 'all'. See EXPERIMENTS for the list.")
    parser.add_argument("--reps-scale", type=float, default=1.0,
                         help="Multiplier on each experiment's default replicate count. "
                              "Use <1 for fast iteration, >1 for publication-grade precision.")
    parser.add_argument("--workers", type=int, default=None, help="Defaults to os.cpu_count().")
    args = parser.parse_args()

    args.iqtree = str(Path(args.iqtree))
    if not os.path.exists(args.iqtree) or not os.access(args.iqtree, os.X_OK):
        raise SystemExit(f"Cannot execute IQ-TREE binary: {args.iqtree}")
    args.outdir.mkdir(parents=True, exist_ok=True)

    selected = list(EXPERIMENTS) if args.which == "all" else [name.strip() for name in args.which.split(",")]
    unknown = [name for name in selected if name not in EXPERIMENTS]
    if unknown:
        raise SystemExit(f"Unknown experiment(s): {unknown}. Choose from: {list(EXPERIMENTS)}")

    print(f"Validating weighted SatuTe formula")
    print(f"  IQ-TREE: {args.iqtree}")
    print(f"  Outdir:  {args.outdir}")
    print(f"  Running: {selected}")

    failures = []
    for name in selected:
        try:
            EXPERIMENTS[name](args)
        except AssertionError as exc:
            failures.append((name, str(exc)))
            print(f"\n  *** {name} FAILED: {exc}")

    section("Summary")
    if failures:
        print(f"{len(failures)}/{len(selected)} experiments failed:")
        for name, message in failures:
            print(f"  - {name}: {message}")
        return 1
    print(f"All {len(selected)} experiments passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
