#!/usr/bin/env python3

import argparse
import math
import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np


STATE_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}
GTR_MODEL = "GTR{1,2,1,1,2,1}+F{0.30,0.20,0.20,0.30}"
NATIVE_REFERENCE_TOL = 2e-4
FORMULA_VARIANTS = ("dominant", "eigenvalue_weighted")


def run(cmd):
    subprocess.run(cmd, check=True)


def parse_fasta(path):
    sequences = {}
    name = None
    chunks = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.upper())
    if name is not None:
        sequences[name] = "".join(chunks)
    return sequences


def parse_iqtree_four_taxon_tree(path):
    text = Path(path).read_text(encoding="utf-8").strip()
    pattern = re.compile(
        r"^\(A:([0-9.eE+-]+),B:([0-9.eE+-]+),"
        r"\(C:([0-9.eE+-]+),D:([0-9.eE+-]+)\):([0-9.eE+-]+)\);$"
    )
    match = pattern.match(text)
    if not match:
        raise ValueError(f"Unsupported generated tree shape in {path}: {text}")
    values = [float(value) for value in match.groups()]
    return {
        "A": values[0],
        "B": values[1],
        "C": values[2],
        "D": values[3],
        "internal": values[4],
    }


def parse_modes(value):
    if not value:
        return []
    return [int(item) for item in value.split(",") if item]


def parse_sat_stat(path):
    header = None
    rows = {}
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            if line.startswith("ID"):
                header = line.rstrip("\n").split("\t")
                continue
            fields = line.rstrip("\n").split("\t")
            if fields and fields[0].isdigit():
                if header is None:
                    raise ValueError(f"No SatuTe header found in {path}")
                row = dict(zip(header, fields))
                if (
                    row.get("Formula") not in FORMULA_VARIANTS
                    or row.get("RateCategory") != "pooled"
                    or row.get("Split") != "A,B"
                ):
                    continue
                rows[row["Formula"]] = {
                    "satC": float(row["satC"]),
                    "satVar": float(row["satVar"]),
                    "satSE": float(row["satSE"]),
                    "satZ": float(row["satZ"]),
                    "satP": float(row["satP"]),
                    "Decision": row["Decision"],
                    "Length": float(row["Length"]),
                    "Modes": parse_modes(row["Modes"]),
                }
    missing = sorted(set(FORMULA_VARIANTS) - set(rows))
    if missing:
        raise ValueError(f"Missing SatuTe pooled rows for {missing} in {path}")
    return rows


def parse_model(model):
    if model == "JC":
        pi = np.array([0.25, 0.25, 0.25, 0.25], dtype=float)
        rates = np.ones(6, dtype=float)
        return rates, pi

    match = re.fullmatch(
        r"GTR\{([^}]+)\}\+F\{([^}]+)\}",
        model,
    )
    if not match:
        raise ValueError(f"Unsupported model syntax: {model}")
    rates = np.array([float(value) for value in match.group(1).split(",")], dtype=float)
    pi = np.array([float(value) for value in match.group(2).split(",")], dtype=float)
    pi /= pi.sum()
    if len(rates) != 6 or len(pi) != 4:
        raise ValueError(f"Unsupported model dimensions: {model}")
    return rates, pi


def build_q(model):
    rates, pi = parse_model(model)
    q = np.zeros((4, 4), dtype=float)
    pairs = [
        (0, 1, rates[0]),
        (0, 2, rates[1]),
        (0, 3, rates[2]),
        (1, 2, rates[3]),
        (1, 3, rates[4]),
        (2, 3, rates[5]),
    ]
    for i, j, rate in pairs:
        q[i, j] = rate * pi[j]
        q[j, i] = rate * pi[i]
    for i in range(4):
        q[i, i] = -np.sum(q[i, :])

    expected_rate = -float(np.dot(pi, np.diag(q)))
    q /= expected_rate
    return q, pi


def reversible_eigendecomposition(q, pi):
    sqrt_pi = np.sqrt(pi)
    inv_sqrt_pi = 1.0 / sqrt_pi
    sym_q = np.diag(sqrt_pi) @ q @ np.diag(inv_sqrt_pi)
    evals, w = np.linalg.eigh(sym_q)
    eigenvectors = np.diag(inv_sqrt_pi) @ w
    inv_eigenvectors = w.T @ np.diag(sqrt_pi)
    return evals, eigenvectors, inv_eigenvectors


def transition_matrix(evals, eigenvectors, inv_eigenvectors, length):
    exp_diag = np.diag(np.exp(evals * length))
    transition = eigenvectors @ exp_diag @ inv_eigenvectors
    transition[transition < 0.0] = np.where(transition[transition < 0.0] > -1e-14, 0.0, transition[transition < 0.0])
    return transition


def mode_indices(evals, formula):
    zero_index = int(np.argmin(np.abs(evals)))
    nonzero = [i for i, value in enumerate(evals) if i != zero_index and abs(value) > 1e-10]
    if formula in {"all_unweighted", "eigenvalue_weighted"}:
        return nonzero
    if formula != "dominant":
        raise ValueError(f"Unknown SatuTe formula: {formula}")

    dominant = max(evals[i] for i in nonzero)
    tol = max(1e-8, abs(dominant) * 1e-6)
    return [i for i in nonzero if abs(evals[i] - dominant) <= tol]


def partial_two_leaf_side(sequences, site, left_leaf, right_leaf, left_length, right_length, transitions):
    partial = np.ones(4, dtype=float)
    for leaf, length in ((left_leaf, left_length), (right_leaf, right_length)):
        symbol = sequences[leaf][site]
        state = STATE_INDEX.get(symbol)
        if state is None:
            continue
        partial *= transitions[length][:, state]
    return partial


def compute_formula(sequences, lengths, model, formula, alpha=0.05):
    q, pi = build_q(model)
    evals, eigenvectors, inv_eigenvectors = reversible_eigendecomposition(q, pi)
    modes = mode_indices(evals, formula)
    if formula == "eigenvalue_weighted":
        dominant = max(evals[modes])
        weights = np.exp((evals[modes] - dominant) * lengths["internal"])
    else:
        weights = np.ones(len(modes), dtype=float)
    transitions = {
        lengths["A"]: transition_matrix(evals, eigenvectors, inv_eigenvectors, lengths["A"]),
        lengths["B"]: transition_matrix(evals, eigenvectors, inv_eigenvectors, lengths["B"]),
        lengths["C"]: transition_matrix(evals, eigenvectors, inv_eigenvectors, lengths["C"]),
        lengths["D"]: transition_matrix(evals, eigenvectors, inv_eigenvectors, lengths["D"]),
    }

    nsites = len(next(iter(sequences.values())))
    coherence = 0.0
    left_second = np.zeros((len(modes), len(modes)), dtype=float)
    right_second = np.zeros((len(modes), len(modes)), dtype=float)
    valid_sites = 0

    for site in range(nsites):
        left_lh = partial_two_leaf_side(
            sequences, site, "A", "B", lengths["A"], lengths["B"], transitions
        )
        right_lh = partial_two_leaf_side(
            sequences, site, "C", "D", lengths["C"], lengths["D"], transitions
        )

        left_post = left_lh * pi
        right_post = right_lh * pi
        left_sum = float(np.sum(left_post))
        right_sum = float(np.sum(right_post))
        if left_sum <= 0.0 or right_sum <= 0.0:
            continue
        left_post /= left_sum
        right_post /= right_sum

        left_factor = np.array(
            [float(np.dot(eigenvectors[:, mode], left_post)) for mode in modes],
            dtype=float,
        )
        right_factor = np.array(
            [float(np.dot(eigenvectors[:, mode], right_post)) for mode in modes],
            dtype=float,
        )

        coherence += float(np.dot(weights, left_factor * right_factor))
        left_second += np.outer(left_factor, left_factor)
        right_second += np.outer(right_factor, right_factor)
        valid_sites += 1

    coherence /= valid_sites
    weight_matrix = np.outer(weights, weights)
    variance = float(np.sum(weight_matrix * (left_second / valid_sites) * (right_second / valid_sites)))
    se = math.sqrt(variance / valid_sites)
    z_score = coherence / se
    p_value = 0.5 * math.erfc(z_score / math.sqrt(2.0))
    return {
        "satC": coherence,
        "satVar": variance,
        "satSE": se,
        "satZ": z_score,
        "satP": p_value,
        "Decision": "informative" if p_value <= alpha else "saturated",
        "Modes": len(modes),
    }


def compute_mode_decomposition(sequences, lengths, model, alpha=0.05):
    q, pi = build_q(model)
    evals, eigenvectors, inv_eigenvectors = reversible_eigendecomposition(q, pi)
    modes = mode_indices(evals, "all_unweighted")
    dominant_modes = set(mode_indices(evals, "dominant"))
    transitions = {
        lengths["A"]: transition_matrix(evals, eigenvectors, inv_eigenvectors, lengths["A"]),
        lengths["B"]: transition_matrix(evals, eigenvectors, inv_eigenvectors, lengths["B"]),
        lengths["C"]: transition_matrix(evals, eigenvectors, inv_eigenvectors, lengths["C"]),
        lengths["D"]: transition_matrix(evals, eigenvectors, inv_eigenvectors, lengths["D"]),
    }

    nsites = len(next(iter(sequences.values())))
    coherence = np.zeros(len(modes), dtype=float)
    left_second = np.zeros(len(modes), dtype=float)
    right_second = np.zeros(len(modes), dtype=float)
    valid_sites = 0

    for site in range(nsites):
        left_lh = partial_two_leaf_side(
            sequences, site, "A", "B", lengths["A"], lengths["B"], transitions
        )
        right_lh = partial_two_leaf_side(
            sequences, site, "C", "D", lengths["C"], lengths["D"], transitions
        )

        left_post = left_lh * pi
        right_post = right_lh * pi
        left_sum = float(np.sum(left_post))
        right_sum = float(np.sum(right_post))
        if left_sum <= 0.0 or right_sum <= 0.0:
            continue
        left_post /= left_sum
        right_post /= right_sum

        for idx, mode in enumerate(modes):
            left_factor = float(np.dot(eigenvectors[:, mode], left_post))
            right_factor = float(np.dot(eigenvectors[:, mode], right_post))
            coherence[idx] += left_factor * right_factor
            left_second[idx] += left_factor * left_factor
            right_second[idx] += right_factor * right_factor
        valid_sites += 1

    rows = []
    for idx, mode in enumerate(modes):
        sat_c = float(coherence[idx] / valid_sites)
        variance = float((left_second[idx] / valid_sites) * (right_second[idx] / valid_sites))
        se = math.sqrt(variance / valid_sites)
        z_score = sat_c / se
        p_value = 0.5 * math.erfc(z_score / math.sqrt(2.0))
        eigenvalue = float(evals[mode])
        half_life = math.log(2.0) / -eigenvalue if eigenvalue < 0.0 else math.inf
        dominant = max(evals[modes])
        decay_weight = math.exp((eigenvalue - dominant) * lengths["internal"])
        rows.append(
            {
                "mode": mode,
                "eigenvalue": eigenvalue,
                "half_life": half_life,
                "decay_weight": decay_weight,
                "satC": sat_c,
                "satSE": se,
                "satZ": z_score,
                "satP": p_value,
                "Decision": "informative" if p_value <= alpha else "saturated",
                "Dominant": mode in dominant_modes,
            }
        )
    return rows


def assert_native_matches_python_reference(label, formula, implementation, reference):
    for field in ("satC", "satVar", "satSE", "satZ", "satP"):
        if not math.isclose(
            implementation[field],
            reference[field],
            rel_tol=0.0,
            abs_tol=NATIVE_REFERENCE_TOL,
        ):
            raise AssertionError(
                f"{label} {formula}: native {field}={implementation[field]} differs from "
                f"Python reference {reference[field]}"
            )
    if implementation["Decision"] != reference["Decision"]:
        raise AssertionError(
            f"{label} {formula}: native decision {implementation['Decision']} differs from "
            f"Python reference {reference['Decision']}"
        )
    if len(implementation["Modes"]) != reference["Modes"]:
        raise AssertionError(
            f"{label} {formula}: native mode count {len(implementation['Modes'])} differs from "
            f"Python reference {reference['Modes']}"
        )


def run_case(iqtree, outdir, label, model, branch_length, seed):
    tree = outdir / f"{label}.tree"
    sim_prefix = outdir / label
    sat_prefix = outdir / f"{label}_sat"
    tree.write_text(f"((A:0.05,B:0.05):{branch_length},C:0.05,D:0.05);\n", encoding="utf-8")

    run([
        iqtree,
        "--alisim",
        str(sim_prefix),
        "-t",
        str(tree),
        "-m",
        model,
        "--length",
        "5000",
        "--seed",
        str(seed),
        "-af",
        "fasta",
        "--quiet",
    ])
    run([
        iqtree,
        "-s",
        str(sim_prefix) + ".fa",
        "-te",
        str(tree),
        "-m",
        model,
        "--satute",
        "--prefix",
        str(sat_prefix),
        "-T",
        "1",
        "--redo",
        "--quiet",
    ])

    sequences = parse_fasta(str(sim_prefix) + ".fa")
    lengths = parse_iqtree_four_taxon_tree(str(sat_prefix) + ".sat.tree")
    implementation = parse_sat_stat(str(sat_prefix) + ".sat.stat")
    references = {
        formula: compute_formula(sequences, lengths, model, formula)
        for formula in FORMULA_VARIANTS
    }
    for formula, reference in references.items():
        assert_native_matches_python_reference(label, formula, implementation[formula], reference)
    mode_decomposition = compute_mode_decomposition(sequences, lengths, model)

    return {
        "label": label,
        "model": model,
        "branch": branch_length,
        "optimized": lengths["internal"],
        "implementation": implementation,
        "references": references,
        "mode_decomposition": mode_decomposition,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Compare SatuTe posterior/eigenvector formula variants on simulated data."
    )
    parser.add_argument("--iqtree", default="build/iqtree3")
    parser.add_argument("--outdir", default="/tmp/iqtree-satute-formula-compare")
    args = parser.parse_args()

    iqtree = str(Path(args.iqtree))
    if not os.path.exists(iqtree) or not os.access(iqtree, os.X_OK):
        raise SystemExit(f"Cannot execute IQ-TREE binary: {iqtree}")

    outdir = Path(args.outdir)
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True)

    cases = [
        ("jc_0.50", "JC", 0.50, 11),
        ("jc_8.00", "JC", 8.00, 800),
        ("gtr_0.50", GTR_MODEL, 0.50, 50),
        ("gtr_4.00", GTR_MODEL, 4.00, 400),
        ("gtr_5.00", GTR_MODEL, 5.00, 500),
        ("gtr_8.00", GTR_MODEL, 8.00, 800),
    ]
    rows = [run_case(iqtree, outdir, *case) for case in cases]

    header = [
        "case",
        "model",
        "sim_branch",
        "opt_branch",
        "dominant_z",
        "dominant_p",
        "dominant_decision",
        "dominant_modes",
        "eigenvalue_weighted_z",
        "eigenvalue_weighted_p",
        "eigenvalue_weighted_decision",
        "eigenvalue_weighted_modes",
    ]
    print("\t".join(header))
    for row in rows:
        dominant = row["implementation"]["dominant"]
        eigenvalue_weighted = row["implementation"]["eigenvalue_weighted"]
        print(
            "\t".join(
                [
                    row["label"],
                    row["model"],
                    f"{row['branch']:.2f}",
                    f"{row['optimized']:.6g}",
                    f"{dominant['satZ']:.6g}",
                    f"{dominant['satP']:.6g}",
                    dominant["Decision"],
                    str(len(dominant["Modes"])),
                    f"{eigenvalue_weighted['satZ']:.6g}",
                    f"{eigenvalue_weighted['satP']:.6g}",
                    eigenvalue_weighted["Decision"],
                    str(len(eigenvalue_weighted["Modes"])),
                ]
            )
        )

    print(
        "\nNative dominant and eigenvalue_weighted rows matched "
        f"the Python references within {NATIVE_REFERENCE_TOL:g}."
    )

    mode_path = outdir / "mode_decomposition.tsv"
    with open(mode_path, "w", encoding="utf-8") as handle:
        handle.write(
            "\t".join(
                [
                    "case",
                    "model",
                    "sim_branch",
                    "opt_branch",
                    "mode",
                    "dominant_mode",
                    "eigenvalue",
                    "half_life",
                    "decay_weight",
                    "mode_satC",
                    "mode_satSE",
                    "mode_z",
                    "mode_p",
                    "mode_decision",
                ]
            )
            + "\n"
        )
        for row in rows:
            for mode_row in row["mode_decomposition"]:
                handle.write(
                    "\t".join(
                        [
                            row["label"],
                            row["model"],
                            f"{row['branch']:.2f}",
                            f"{row['optimized']:.10g}",
                            str(mode_row["mode"]),
                            "yes" if mode_row["Dominant"] else "no",
                            f"{mode_row['eigenvalue']:.10g}",
                            f"{mode_row['half_life']:.10g}",
                            f"{mode_row['decay_weight']:.10g}",
                            f"{mode_row['satC']:.10g}",
                            f"{mode_row['satSE']:.10g}",
                            f"{mode_row['satZ']:.10g}",
                            f"{mode_row['satP']:.10g}",
                            mode_row["Decision"],
                        ]
                    )
                    + "\n"
                )

    print(f"\nMode decomposition written to {mode_path}")
    print(f"Outputs are in {outdir}")


if __name__ == "__main__":
    main()
