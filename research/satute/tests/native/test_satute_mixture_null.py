#!/usr/bin/env python3
"""Check native mixture scores against directly built centred mixture scores."""

import argparse
import csv
import itertools
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT / "src"))
from satute_analysis.models import build_q, reversible_eigendecomposition

MODEL = "GTR{1,2,1,1,2,1}+F{0.30,0.20,0.20,0.30}"
RATES = np.array([0.03338775338, 0.2519159176, 0.8202684819, 2.894427847])
PATTERNS = tuple(itertools.product(range(4), repeat=4))


def table(path):
    return list(csv.DictReader(
        (line for line in path.read_text().splitlines() if line and not line.startswith("#")),
        delimiter="\t",
    ))


def run_native(iqtree, alignment, tree, model, prefix, *extra, disable_cache=False, threads=1):
    environment = dict(os.environ)
    if disable_cache:
        environment["SATUTE_DISABLE_DIRECTIONAL_CACHE"] = "1"
    else:
        environment.pop("SATUTE_DISABLE_DIRECTIONAL_CACHE", None)
    subprocess.run([
        str(iqtree), "-s", str(alignment), "-te", str(tree), "-m", model,
        "-blfix", "--satute", "--satute-sites", "--prefix", str(prefix),
        "-T", str(threads), "--redo", "--quiet", *extra,
    ], check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=environment)
    rows = table(Path(str(prefix) + ".sat.stat"))
    pooled = {r["Formula"]: r for r in rows if r["Split"] == "A,B" and r["RateCategory"] == "pooled"}
    assert set(pooled) == {"dominant", "eigenvalue_weighted"}, pooled
    return pooled, rows, table(Path(str(prefix) + ".sat.site"))


def direct_score(pinv, length, model):
    """Centred mixture score per pattern, built from explicit transition
    matrices instead of the native recursion.

    Each positive-rate category contributes the centred mode product
    sum_i w_i (L_i - mean L_i)(R_i - mean R_i), weighted by that category's
    responsibility. Centring is what keeps the statistic a covariance when
    the fitted model does not generate the data -- here the alignment is
    uniform over all 256 patterns, which no GTR model produces.
    """
    q, pi = build_q(model)
    evals, eigenvectors, _ = reversible_eigendecomposition(q, pi)
    zero_index = int(np.argmin(np.abs(evals)))
    modes = [i for i in range(len(evals)) if i != zero_index and abs(evals[i]) > 1e-10]
    dominant = max(evals[modes])
    rates = RATES / (1.0 - pinv)
    shift = dominant * min(rates) * length
    null = np.zeros(len(PATTERNS))
    weighted = np.zeros(len(PATTERNS))
    for rate in ([0.0] if pinv else []) + list(rates):
        proportion = pinv if rate == 0.0 else (1.0 - pinv) / 4.0
        side = expm(q * rate * 0.05)
        left = np.array([side[:, p[0]] * side[:, p[1]] for p in PATTERNS])
        right = np.array([side[:, p[2]] * side[:, p[3]] for p in PATTERNS])
        left_mass, right_mass = left @ pi, right @ pi
        if rate == 0.0:
            # No branch signal, but the joint probability stays in the null.
            density = proportion * (left * pi * right).sum(axis=1)
            score = np.zeros(len(PATTERNS))
        else:
            density = proportion * left_mass * right_mass
            left_factor = (left * pi / left_mass[:, None]) @ eigenvectors[:, modes]
            right_factor = (right * pi / right_mass[:, None]) @ eigenvectors[:, modes]
            left_factor -= left_factor.mean(axis=0)
            right_factor -= right_factor.mean(axis=0)
            score = (left_factor * right_factor) @ np.exp(evals[modes] * rate * length - shift)
        null += density
        weighted += density * score
    assert np.isclose(null.sum(), 1.0)
    return weighted / null


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("iqtree", type=Path)
    parser.add_argument("outdir", type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    alignment = args.outdir / "patterns.fa"
    alignment.write_text("".join(
        f">{taxon}\n" + "".join("ACGT"[p[i]] for p in PATTERNS) + "\n"
        for i, taxon in enumerate("ABCD")
    ))
    tree = args.outdir / "tree.nwk"
    tree.write_text("(A:0.05,B:0.05,(C:0.05,D:0.05):3);\n")
    failures = []
    for base, pinv in ((MODEL, 0.0), (MODEL, 0.1), ("JC", 0.0), ("JC", 0.1)):
        name = ("jc_" if base == "JC" else "gtr_") + ("gamma" if pinv == 0.0 else "invar_gamma")
        model = base + ("+I{0.1}" if pinv else "") + "+G4{0.5}"
        pooled, rows, sites = run_native(args.iqtree, alignment, tree, model, args.outdir / name)
        row = pooled["eigenvalue_weighted"]
        actual = np.array([float(r["SiteCoherence"]) for r in sites
                           if r["ID"] == row["ID"] and r["Formula"] == "eigenvalue_weighted"])
        expected = direct_score(pinv, 3.0, base)
        assert len(actual) == len(expected)
        error = float(np.max(np.abs(actual - expected)))
        if error > 2e-7:
            failures.append(f"{name}: native site score differs from the centred mixture score by {error:g}")
        expected_variance = float(np.var(expected, ddof=1))
        if not math.isclose(float(row["satVar"]), expected_variance, rel_tol=2e-7, abs_tol=2e-7):
            failures.append(f"{name}: pooled variance does not studentize the mixture site scores")
        if not math.isclose(float(row["satC"]), float(expected.mean()), rel_tol=2e-7, abs_tol=2e-7):
            failures.append(f"{name}: pooled mean does not match the centred mixture score")
        for diagnostic in rows:
            if diagnostic["RateCategory"] != "pooled" and diagnostic["satP"] != "NA":
                failures.append(f"{name}: hard category {diagnostic['RateCategory']} reports a nominal p-value")
                break
        print(f"{name}: centred mixture site-score maximum error={error:g}")
        if base == "JC":
            for field in ("satC", "satVar", "satZ"):
                assert math.isclose(float(pooled["dominant"][field]), float(row[field]), rel_tol=2e-10, abs_tol=2e-10)
        if base == MODEL and pinv == 0.1:
            for label, options in (("fallback", {"disable_cache": True}), ("threads", {"threads": 2})):
                alternative, _, exported = run_native(args.iqtree, alignment, tree, model, args.outdir / label, **options)
                for formula in pooled:
                    for field in ("satC", "satVar", "satZ"):
                        assert math.isclose(float(pooled[formula][field]), float(alternative[formula][field]), rel_tol=2e-10, abs_tol=2e-10), (label, formula, field)
                compared = [float(r["SiteCoherence"]) for r in exported
                            if r["ID"] == alternative["eigenvalue_weighted"]["ID"] and r["Formula"] == "eigenvalue_weighted"]
                np.testing.assert_allclose(actual, compared, rtol=2e-10, atol=2e-10)

    signal = args.outdir / "signal.fa"
    signal_patterns = PATTERNS + tuple((i, i, i, i) for i in range(4)) * 16
    signal.write_text("".join(f">{taxon}\n" + "".join("ACGT"[p[i]] for p in signal_patterns) + "\n"
                              for i, taxon in enumerate("ABCD")))
    plain, _, _ = run_native(args.iqtree, signal, tree, MODEL, args.outdir / "plain")
    enhanced, _, _ = run_native(args.iqtree, signal, tree, MODEL, args.outdir / "enhanced", "--satute-enhanced-weights")
    for field in ("satC", "satVar", "satZ"):
        if not math.isclose(float(plain["eigenvalue_weighted"][field]), float(enhanced["eigenvalue_weighted"][field]), rel_tol=2e-7, abs_tol=2e-7):
            failures.append(f"enhanced: full-rank model-mean projection changed {field}")
    if failures:
        raise AssertionError("\n".join(failures))
    print("Native mixture-null and model-mean weighting regression passed")


if __name__ == "__main__":
    main()
