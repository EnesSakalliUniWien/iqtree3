#!/usr/bin/env python3

import argparse
import math
import os
import re
import shutil
import subprocess
from collections import Counter
from pathlib import Path


STATE_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}
JC_FREQUENCIES = [0.25, 0.25, 0.25, 0.25]
JC_EIGENVECTORS = [
    [math.sqrt(2.0), math.sqrt(2.0 / 3.0), math.sqrt(1.0 / 3.0)],
    [-math.sqrt(2.0), math.sqrt(2.0 / 3.0), math.sqrt(1.0 / 3.0)],
    [0.0, -2.0 * math.sqrt(2.0 / 3.0), math.sqrt(1.0 / 3.0)],
    [0.0, 0.0, -math.sqrt(3.0)],
]
NUMERIC_TOL = 2e-4
SIM_MODEL = "JC+G4{0.5}"
GTR_GAMMA_MODEL = "GTR{1,2,1,1,2,1}+F{0.30,0.20,0.20,0.30}+G4{0.5}"
TREE = "((A:0.10,B:0.10):0.40,C:0.10,D:0.10);\n"

SUPPORTED_MODELS = [
    ("gamma", "JC+G4{0.5}"),
    ("freerate", "JC+R4"),
    ("invar_gamma", "JC+I{0.1}+G4{0.5}"),
    ("invar_freerate", "JC+I+R4"),
    ("gtr_gamma", GTR_GAMMA_MODEL),
]

PURE_INVAR_MODEL = "JC+I{0.1}"
PURE_INVAR_ERROR = "IQ-TREE does not emit per-site category rows for a pure +I model"


def run(cmd, expect_success=True):
    completed = subprocess.run(
        cmd,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if expect_success and completed.returncode != 0:
        raise RuntimeError(
            "Command failed with exit code "
            f"{completed.returncode}:\n{' '.join(cmd)}\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
    if not expect_success and completed.returncode == 0:
        raise RuntimeError(f"Command unexpectedly succeeded: {' '.join(cmd)}")
    return completed


def parse_rate_file(path):
    header = None
    counts = Counter()
    category_rates = {}
    site_categories = {}

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[0] == "Site":
                header = fields
                continue
            if header is None:
                raise ValueError(f"No .rate header found in {path}")
            row = dict(zip(header, fields))
            cat = row["Cat"]
            counts[cat] += 1
            category_rates.setdefault(cat, float(row["C_Rate"]))
            site_categories[int(row["Site"]) - 1] = cat

    return counts, category_rates, site_categories


def parse_sat_stat(path, split=None, formula="dominant"):
    header = None
    categories = {}
    pooled_row = None

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[0] == "ID":
                header = fields
                continue
            if header is None:
                raise ValueError(f"No SatuTe header found in {path}")
            if not fields[0].isdigit():
                continue
            row = dict(zip(header, fields))
            if row["Formula"] != formula:
                continue
            if split is not None and row["Split"] != split:
                continue
            if row["RateCategory"] == "pooled":
                pooled_row = row
                continue
            if row.get("FDR_BY") != "NA" or row.get("DecisionFDR") != "not_tested":
                raise AssertionError(
                    f"Rate-category row entered the pooled {formula} FDR family in {path}"
                )
            cat = row["RateCategory"]
            categories.setdefault(
                cat,
                {
                    "sites": int(row["RateSites"]),
                    "rate": float(row["RateMultiplier"]),
                    "length": float(row["Length"]),
                    "effective_length": float(row["EffectiveLength"]),
                    "satC": float_or_none(row["satC"]),
                    "satVar": float_or_none(row["satVar"]),
                    "satSE": float_or_none(row["satSE"]),
                    "satZ": float_or_none(row["satZ"]),
                    "valid_sites": int(row["ValidSites"]),
                },
            )

    if pooled_row is None:
        raise ValueError(f"No {formula} pooled row found in {path}")
    if not categories:
        raise ValueError(f"No {formula} category rows found in {path}")
    return categories, pooled_row


def float_or_none(value):
    if value == "NA":
        return None
    return float(value)


def parse_fasta(path):
    sequences = {}
    name = None
    chunks = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
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


def parse_four_taxon_sat_tree(path):
    text = path.read_text(encoding="utf-8").strip()
    pattern = re.compile(
        r"^\(A:([0-9.eE+-]+),B:([0-9.eE+-]+),"
        r"\(C:([0-9.eE+-]+),D:([0-9.eE+-]+)\):([0-9.eE+-]+)\);$"
    )
    match = pattern.match(text)
    if not match:
        raise ValueError(f"Unsupported four-taxon SatuTe tree shape in {path}: {text}")
    values = [float(value) for value in match.groups()]
    return {
        "A": values[0],
        "B": values[1],
        "C": values[2],
        "D": values[3],
        "internal": values[4],
    }


def jc_transition(length):
    same = 0.25 + 0.75 * math.exp(-4.0 * length / 3.0)
    diff = 0.25 - 0.25 * math.exp(-4.0 * length / 3.0)
    return [
        [same if row == col else diff for col in range(4)]
        for row in range(4)
    ]


def parse_gtr_components(model):
    match = re.search(r"GTR\{([^}]+)\}\+F\{([^}]+)\}", model)
    if not match:
        raise ValueError(f"Unsupported GTR model syntax: {model}")
    rates = [float(value) for value in match.group(1).split(",")]
    frequencies = [float(value) for value in match.group(2).split(",")]
    total = sum(frequencies)
    frequencies = [value / total for value in frequencies]
    if len(rates) != 6 or len(frequencies) != 4:
        raise ValueError(f"Unsupported GTR model dimensions: {model}")
    return rates, frequencies


def build_gtr_q(model):
    rates, frequencies = parse_gtr_components(model)
    q = [[0.0 for _ in range(4)] for _ in range(4)]
    pairs = [
        (0, 1, rates[0]),
        (0, 2, rates[1]),
        (0, 3, rates[2]),
        (1, 2, rates[3]),
        (1, 3, rates[4]),
        (2, 3, rates[5]),
    ]
    for i, j, rate in pairs:
        q[i][j] = rate * frequencies[j]
        q[j][i] = rate * frequencies[i]
    for i in range(4):
        q[i][i] = -sum(q[i])
    expected_rate = -sum(frequencies[i] * q[i][i] for i in range(4))
    return [[value / expected_rate for value in row] for row in q], frequencies


def jacobi_eigen_symmetric(matrix):
    n = len(matrix)
    a = [row[:] for row in matrix]
    vectors = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

    for _iteration in range(100):
        p, q = 0, 1
        max_offdiag = abs(a[p][q])
        for i in range(n):
            for j in range(i + 1, n):
                value = abs(a[i][j])
                if value > max_offdiag:
                    max_offdiag = value
                    p, q = i, j
        if max_offdiag < 1e-14:
            break

        if abs(a[p][p] - a[q][q]) < 1e-30:
            angle = math.pi / 4.0
        else:
            angle = 0.5 * math.atan2(2.0 * a[p][q], a[q][q] - a[p][p])
        c = math.cos(angle)
        s = math.sin(angle)

        app = c * c * a[p][p] - 2.0 * s * c * a[p][q] + s * s * a[q][q]
        aqq = s * s * a[p][p] + 2.0 * s * c * a[p][q] + c * c * a[q][q]
        a[p][p] = app
        a[q][q] = aqq
        a[p][q] = 0.0
        a[q][p] = 0.0

        for k in range(n):
            if k == p or k == q:
                continue
            akp = c * a[k][p] - s * a[k][q]
            akq = s * a[k][p] + c * a[k][q]
            a[k][p] = akp
            a[p][k] = akp
            a[k][q] = akq
            a[q][k] = akq

        for k in range(n):
            vkp = c * vectors[k][p] - s * vectors[k][q]
            vkq = s * vectors[k][p] + c * vectors[k][q]
            vectors[k][p] = vkp
            vectors[k][q] = vkq

    order = sorted(range(n), key=lambda i: a[i][i])
    eigenvalues = [a[i][i] for i in order]
    eigenvectors = [[vectors[row][i] for i in order] for row in range(n)]
    return eigenvalues, eigenvectors


def gtr_eigendecomposition(model):
    q, frequencies = build_gtr_q(model)
    sqrt_pi = [math.sqrt(value) for value in frequencies]
    inv_sqrt_pi = [1.0 / value for value in sqrt_pi]
    sym_q = [
        [sqrt_pi[i] * q[i][j] * inv_sqrt_pi[j] for j in range(4)]
        for i in range(4)
    ]
    eigenvalues, sym_eigenvectors = jacobi_eigen_symmetric(sym_q)
    eigenvectors = [
        [inv_sqrt_pi[state] * sym_eigenvectors[state][mode] for mode in range(4)]
        for state in range(4)
    ]
    inv_eigenvectors = [
        [sym_eigenvectors[mode_state][mode] * sqrt_pi[mode_state] for mode_state in range(4)]
        for mode in range(4)
    ]
    return eigenvalues, eigenvectors, inv_eigenvectors, frequencies


def gtr_transition(eigenvalues, eigenvectors, inv_eigenvectors, length):
    transition = [[0.0 for _ in range(4)] for _ in range(4)]
    for parent_state in range(4):
        for child_state in range(4):
            value = 0.0
            for mode in range(4):
                value += (
                    eigenvectors[parent_state][mode]
                    * math.exp(eigenvalues[mode] * length)
                    * inv_eigenvectors[mode][child_state]
                )
            transition[parent_state][child_state] = value
    return transition


def dominant_modes(eigenvalues):
    zero_index = min(range(len(eigenvalues)), key=lambda i: abs(eigenvalues[i]))
    nonzero = [
        i for i, value in enumerate(eigenvalues)
        if i != zero_index and abs(value) > 1e-10
    ]
    largest = max(eigenvalues[i] for i in nonzero)
    tol = max(1e-8, abs(largest) * 1e-6)
    return [i for i in nonzero if abs(eigenvalues[i] - largest) <= tol]


def all_nonzero_modes(eigenvalues):
    zero_index = min(range(len(eigenvalues)), key=lambda i: abs(eigenvalues[i]))
    return [
        i for i, value in enumerate(eigenvalues)
        if i != zero_index and abs(value) > 1e-10
    ]


def formula_modes_and_weights(eigenvalues, formula, effective_branch_length):
    if formula == "dominant":
        modes = dominant_modes(eigenvalues)
        weights = [1.0 for _mode in modes]
    elif formula == "all_unweighted":
        modes = all_nonzero_modes(eigenvalues)
        weights = [1.0 for _mode in modes]
    elif formula == "eigenvalue_weighted":
        modes = all_nonzero_modes(eigenvalues)
        dominant = max(eigenvalues[mode] for mode in modes)
        weights = [
            math.exp((eigenvalues[mode] - dominant) * effective_branch_length)
            for mode in modes
        ]
    else:
        raise ValueError(f"Unsupported SatuTe formula: {formula}")
    return modes, weights


def side_partial(sequences, site, leaves, lengths, rate_multiplier, transition_builder):
    partial = [1.0, 1.0, 1.0, 1.0]
    for leaf in leaves:
        state = STATE_INDEX.get(sequences[leaf][site])
        if state is None:
            continue
        transition = transition_builder(lengths[leaf] * rate_multiplier)
        for parent_state in range(4):
            partial[parent_state] *= transition[parent_state][state]
    return partial


def posterior_from_partial(partial, frequencies):
    posterior = [partial[state] * frequencies[state] for state in range(4)]
    total = sum(posterior)
    if not math.isfinite(total) or abs(total) <= 0.0:
        return None
    return [value / total for value in posterior]


def mode_factors(posterior, eigenvectors, modes):
    return [
        sum(eigenvectors[state][mode] * posterior[state] for state in range(4))
        for mode in modes
    ]


def compute_reference(
    sequences,
    lengths,
    sites,
    rate_multiplier,
    frequencies,
    eigenvectors,
    modes,
    mode_weights,
    transition_builder,
):
    coherence_sum = 0.0
    nmodes = len(modes)
    left_second = [[0.0 for _ in range(nmodes)] for _ in range(nmodes)]
    right_second = [[0.0 for _ in range(nmodes)] for _ in range(nmodes)]
    valid_sites = 0

    for site in sites:
        left_post = posterior_from_partial(
            side_partial(sequences, site, ("A", "B"), lengths, rate_multiplier, transition_builder),
            frequencies,
        )
        right_post = posterior_from_partial(
            side_partial(sequences, site, ("C", "D"), lengths, rate_multiplier, transition_builder),
            frequencies,
        )
        if left_post is None or right_post is None:
            continue

        left = mode_factors(left_post, eigenvectors, modes)
        right = mode_factors(right_post, eigenvectors, modes)
        coherence_sum += sum(mode_weights[idx] * left[idx] * right[idx] for idx in range(nmodes))
        for i in range(nmodes):
            for j in range(nmodes):
                left_second[i][j] += left[i] * left[j]
                right_second[i][j] += right[i] * right[j]
        valid_sites += 1

    if valid_sites == 0:
        raise AssertionError("No valid sites for Python JC reference calculation")

    sat_c = coherence_sum / valid_sites
    variance = 0.0
    for i in range(nmodes):
        for j in range(nmodes):
            variance += (
                mode_weights[i]
                * mode_weights[j]
                * (left_second[i][j] / valid_sites)
                * (right_second[i][j] / valid_sites)
            )
    sat_se = math.sqrt(variance / valid_sites)
    sat_z = sat_c / sat_se
    return {
        "satC": sat_c,
        "satVar": variance,
        "satSE": sat_se,
        "satZ": sat_z,
        "valid_sites": valid_sites,
    }


def compute_jc_reference(sequences, lengths, sites, rate_multiplier):
    return compute_reference(
        sequences,
        lengths,
        sites,
        rate_multiplier,
        JC_FREQUENCIES,
        JC_EIGENVECTORS,
        [0, 1, 2],
        [1.0, 1.0, 1.0],
        jc_transition,
    )


def compute_gtr_reference(sequences, lengths, sites, rate_multiplier, model, formula):
    eigenvalues, eigenvectors, inv_eigenvectors, frequencies = gtr_eigendecomposition(model)
    modes, weights = formula_modes_and_weights(
        eigenvalues,
        formula,
        lengths["internal"] * rate_multiplier,
    )
    return compute_reference(
        sequences,
        lengths,
        sites,
        rate_multiplier,
        frequencies,
        eigenvectors,
        modes,
        weights,
        lambda length: gtr_transition(eigenvalues, eigenvectors, inv_eigenvectors, length),
    )


def assert_close(label, field, observed, expected):
    if observed is None:
        raise AssertionError(f"{label}: native {field} is NA but Python reference is {expected}")
    if not math.isclose(observed, expected, rel_tol=0.0, abs_tol=NUMERIC_TOL):
        raise AssertionError(
            f"{label}: native {field}={observed} differs from Python reference {expected}"
        )


def assert_category_consistency(label, rate_file, sat_file):
    rate_counts, rate_values, _site_categories = parse_rate_file(rate_file)
    sat_categories, _pooled_row = parse_sat_stat(sat_file)

    if not rate_counts:
        raise AssertionError(f"{label}: IQ-TREE .rate file contains no site rows")

    for cat, count in rate_counts.items():
        if cat not in sat_categories:
            raise AssertionError(f"{label}: SatuTe is missing IQ-TREE category {cat}")
        sat = sat_categories[cat]
        if sat["sites"] != count:
            raise AssertionError(
                f"{label}: category {cat} has {sat['sites']} SatuTe sites, "
                f"but IQ-TREE .rate has {count}"
            )
        if not math.isclose(sat["rate"], rate_values[cat], rel_tol=0.0, abs_tol=1e-4):
            raise AssertionError(
                f"{label}: category {cat} rate {sat['rate']} does not match "
                f"IQ-TREE .rate {rate_values[cat]}"
            )
        expected_effective = sat["length"] * sat["rate"]
        if not math.isclose(
            sat["effective_length"], expected_effective, rel_tol=0.0, abs_tol=1e-8
        ):
            raise AssertionError(
                f"{label}: category {cat} effective length {sat['effective_length']} "
                f"does not equal length*rate {expected_effective}"
            )

    extras = {
        cat: row["sites"]
        for cat, row in sat_categories.items()
        if cat not in rate_counts
    }
    nonzero_extras = {cat: sites for cat, sites in extras.items() if sites != 0}
    if nonzero_extras:
        raise AssertionError(
            f"{label}: SatuTe has nonzero category rows absent from IQ-TREE .rate: "
            f"{nonzero_extras}"
        )

    return rate_counts, extras


def assert_jc_gamma_numeric_reference(alignment, prefix):
    sequences = parse_fasta(alignment)
    lengths = parse_four_taxon_sat_tree(prefix.with_suffix(".sat.tree"))
    rate_counts, rate_values, site_categories = parse_rate_file(prefix.with_suffix(".rate"))
    sat_categories, pooled_row = parse_sat_stat(prefix.with_suffix(".sat.stat"), split="A,B")

    category_refs = {}
    valid_total = 0
    pooled_c = 0.0
    pooled_var = 0.0

    for cat, count in rate_counts.items():
        sites = [site for site, site_cat in site_categories.items() if site_cat == cat]
        if len(sites) != count:
            raise AssertionError(f"gamma: category {cat} site list/count mismatch")
        reference = compute_jc_reference(sequences, lengths, sites, rate_values[cat])
        category_refs[cat] = reference
        valid_total += reference["valid_sites"]

        native = sat_categories[cat]
        if native["valid_sites"] != reference["valid_sites"]:
            raise AssertionError(
                f"gamma category {cat}: native valid sites {native['valid_sites']} "
                f"differs from Python reference {reference['valid_sites']}"
            )
        for field in ("satC", "satVar", "satSE", "satZ"):
            assert_close(f"gamma category {cat}", field, native[field], reference[field])

    if valid_total != sum(rate_counts.values()):
        raise AssertionError("gamma: pooled reference site count differs from IQ-TREE .rate count")

    for reference in category_refs.values():
        weight = reference["valid_sites"] / valid_total
        pooled_c += weight * reference["satC"]
        pooled_var += weight * reference["satVar"]
    pooled_se = math.sqrt(pooled_var / valid_total)
    pooled_z = pooled_c / pooled_se

    pooled_expected = {
        "satC": pooled_c,
        "satVar": pooled_var,
        "satSE": pooled_se,
        "satZ": pooled_z,
    }
    for field, expected in pooled_expected.items():
        assert_close("gamma pooled", field, float_or_none(pooled_row[field]), expected)

    return sorted(category_refs)


def assert_gtr_gamma_numeric_reference(alignment, prefix, model):
    sequences = parse_fasta(alignment)
    lengths = parse_four_taxon_sat_tree(prefix.with_suffix(".sat.tree"))
    rate_counts, _rate_values, site_categories = parse_rate_file(prefix.with_suffix(".rate"))
    checked = {}

    for formula in ("dominant", "eigenvalue_weighted"):
        sat_categories, pooled_row = parse_sat_stat(
            prefix.with_suffix(".sat.stat"),
            split="A,B",
            formula=formula,
        )
        category_refs = {}
        valid_total = 0
        pooled_c = 0.0
        pooled_var = 0.0

        for cat, count in rate_counts.items():
            sites = [site for site, site_cat in site_categories.items() if site_cat == cat]
            if len(sites) != count:
                raise AssertionError(f"gtr_gamma {formula}: category {cat} site list/count mismatch")
            reference = compute_gtr_reference(
                sequences,
                lengths,
                sites,
                sat_categories[cat]["rate"],
                model,
                formula,
            )
            category_refs[cat] = reference
            valid_total += reference["valid_sites"]

            native = sat_categories[cat]
            if native["valid_sites"] != reference["valid_sites"]:
                raise AssertionError(
                    f"gtr_gamma {formula} category {cat}: native valid sites {native['valid_sites']} "
                    f"differs from Python reference {reference['valid_sites']}"
                )
            for field in ("satC", "satVar", "satSE", "satZ"):
                assert_close(f"gtr_gamma {formula} category {cat}", field, native[field], reference[field])

        if valid_total != sum(rate_counts.values()):
            raise AssertionError(f"gtr_gamma {formula}: pooled reference site count differs from IQ-TREE .rate count")

        for reference in category_refs.values():
            weight = reference["valid_sites"] / valid_total
            pooled_c += weight * reference["satC"]
            pooled_var += weight * reference["satVar"]
        pooled_se = math.sqrt(pooled_var / valid_total)
        pooled_z = pooled_c / pooled_se

        pooled_expected = {
            "satC": pooled_c,
            "satVar": pooled_var,
            "satSE": pooled_se,
            "satZ": pooled_z,
        }
        for field, expected in pooled_expected.items():
            assert_close(f"gtr_gamma {formula} pooled", field, float_or_none(pooled_row[field]), expected)

        checked[formula] = sorted(category_refs)

    return checked


def main():
    parser = argparse.ArgumentParser(
        description="Verify native SatuTe rate categories against IQ-TREE .rate output."
    )
    parser.add_argument("iqtree", nargs="?", default="build/iqtree3")
    parser.add_argument("outdir", nargs="?", default="/tmp/iqtree-satute-rate-categories")
    args = parser.parse_args()

    iqtree = Path(args.iqtree)
    if not iqtree.is_file():
        raise SystemExit(f"Cannot find IQ-TREE binary: {iqtree}")
    if not os.access(iqtree, os.X_OK):
        raise SystemExit(f"Cannot execute IQ-TREE binary: {iqtree}")

    outdir = Path(args.outdir)
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True)

    tree = outdir / "rate_test.tree"
    tree.write_text(TREE, encoding="utf-8")

    sim_prefix = outdir / "sim"
    run(
        [
            str(iqtree),
            "--alisim",
            str(sim_prefix),
            "-t",
            str(tree),
            "-m",
            SIM_MODEL,
            "--length",
            "1000",
            "--seed",
            "20260619",
            "-af",
            "fasta",
            "--quiet",
        ]
    )

    alignment = sim_prefix.with_suffix(".fa")
    print("case\tmodel\tmatched_categories\textra_zero_categories")
    for label, model in SUPPORTED_MODELS:
        case_dir = outdir / label
        case_dir.mkdir()
        prefix = case_dir / "example"
        run(
            [
                str(iqtree),
                "-s",
                str(alignment),
                "-te",
                str(tree),
                "-m",
                model,
                "--prefix",
                str(prefix),
                "-T",
                "1",
                "--redo",
                "--quiet",
                "--rate",
                "--satute",
            ]
        )
        rate_counts, extras = assert_category_consistency(
            label, prefix.with_suffix(".rate"), prefix.with_suffix(".sat.stat")
        )
        numeric_checked = ""
        if label == "gamma":
            numeric_categories = assert_jc_gamma_numeric_reference(alignment, prefix)
            numeric_checked = f"; numeric_reference={numeric_categories}"
        if label == "gtr_gamma":
            numeric_categories = assert_gtr_gamma_numeric_reference(alignment, prefix, model)
            numeric_checked = f"; numeric_reference={numeric_categories}"
        extra_zero = sorted(cat for cat, sites in extras.items() if sites == 0)
        print(
            f"{label}\t{model}\t{dict(sorted(rate_counts.items(), key=lambda item: int(item[0])))}"
            f"\t{extra_zero}{numeric_checked}"
        )

    pure_dir = outdir / "pure_invar"
    pure_dir.mkdir()
    pure_prefix = pure_dir / "example"
    failed = run(
        [
            str(iqtree),
            "-s",
            str(alignment),
            "-te",
            str(tree),
            "-m",
            PURE_INVAR_MODEL,
            "--prefix",
            str(pure_prefix),
            "-T",
            "1",
            "--redo",
            "--quiet",
            "--rate",
            "--satute",
        ],
        expect_success=False,
    )
    combined_output = failed.stdout + failed.stderr
    if PURE_INVAR_ERROR not in combined_output:
        raise AssertionError(
            "Pure +I failed, but not with the expected IQ-TREE category message:\n"
            + combined_output
        )
    print(f"pure_invar\t{PURE_INVAR_MODEL}\texpected_failure\tpure +I has no IQ-TREE category rows")
    print(f"SatuTe rate-category regression passed; outputs are in {outdir}")


if __name__ == "__main__":
    main()
