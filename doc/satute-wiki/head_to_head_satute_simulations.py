#!/usr/bin/env python3

import argparse
import csv
import math
import os
import random
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

try:
    import numpy as np
except ModuleNotFoundError as exc:
    raise SystemExit(
        "This script requires NumPy. Run it with a Python environment that has NumPy installed, "
        "or set PYTHON_BIN to that interpreter before launching the documented commands."
    ) from exc


STATE_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}
PAPER_BRANCH_LENGTHS = "0.1,0.2,0.3,0.4,0.5,0.8,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,7.5,10.0"
GTR_PF06346_MODEL = "GTR{0.6676,3.7807,4.2833,0.5354,0.8718,1.0}+F{0.125,0.436,0.191,0.245}"
MODEL_ALIASES = {
    "JC": "JC",
    "GTR_PF06346": GTR_PF06346_MODEL,
    "GTR_EvoNAPS_PF06346": GTR_PF06346_MODEL,
}


class Node:
    next_id = 0

    def __init__(self, name=""):
        self.id = Node.next_id
        Node.next_id += 1
        self.name = name
        self.neighbors = []

    def add_neighbor(self, other, length):
        self.neighbors.append((other, length))


def run(cmd):
    subprocess.run(cmd, check=True)


def parse_csv_numbers(text, cast=float):
    return [cast(value) for value in text.split(",") if value.strip()]


def resolve_model_alias(alias):
    if alias in MODEL_ALIASES:
        return MODEL_ALIASES[alias]
    return alias


def sniff_delimiter(path):
    first = Path(path).read_text(encoding="utf-8").splitlines()[0]
    return "\t" if first.count("\t") >= first.count(",") else ","


def load_branch_length_pools(path):
    pools = {"internal": [], "external": []}
    if not path:
        return pools

    delimiter = sniff_delimiter(path)
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        lowered = {name.lower(): name for name in (reader.fieldnames or [])}
        length_col = lowered.get("length") or lowered.get("branch_length") or lowered.get("blen") or lowered.get("bl")
        kind_col = lowered.get("kind") or lowered.get("type") or lowered.get("branch_type")
        if not length_col or not kind_col:
            raise ValueError(
                "EvoNAPS branch-length table must have a length column "
                "(length, branch_length, blen, or bl) and a type column "
                "(kind, type, or branch_type)."
            )
        for row in reader:
            kind = row[kind_col].strip().lower()
            if kind == "i" or kind.startswith("int"):
                kind = "internal"
            elif kind == "e" or kind.startswith("ext") or kind.startswith("term"):
                kind = "external"
            else:
                continue
            length = float(row[length_col])
            if 0.04 <= length <= 0.2:
                pools[kind].append(length)

    if not pools["internal"] or not pools["external"]:
        raise ValueError(f"No usable internal/external branch lengths found in {path}")
    return pools


def sample_length(rng, pools, kind, fallback):
    if pools.get(kind):
        return rng.choice(pools[kind])
    return fallback


def add_edge(a, b, length):
    a.add_neighbor(b, length)
    b.add_neighbor(a, length)


def parse_newick(path):
    text = Path(path).read_text(encoding="utf-8").strip()
    if text.endswith(";"):
        text = text[:-1]
    Node.next_id = 0
    idx = 0

    def skip_ws():
        nonlocal idx
        while idx < len(text) and text[idx].isspace():
            idx += 1

    def parse_label():
        nonlocal idx
        skip_ws()
        start = idx
        while idx < len(text) and text[idx] not in ":,()":
            idx += 1
        return text[start:idx].strip()

    def parse_length():
        nonlocal idx
        skip_ws()
        if idx >= len(text) or text[idx] != ":":
            return 0.0
        idx += 1
        skip_ws()
        start = idx
        while idx < len(text) and text[idx] not in ",()":
            idx += 1
        return float(text[start:idx])

    def parse_subtree():
        nonlocal idx
        skip_ws()
        if text[idx] == "(":
            idx += 1
            node = Node()
            while True:
                child, child_length = parse_subtree()
                add_edge(node, child, child_length)
                skip_ws()
                if idx < len(text) and text[idx] == ",":
                    idx += 1
                    continue
                if idx < len(text) and text[idx] == ")":
                    idx += 1
                    break
                raise ValueError(f"Unexpected Newick character near {text[idx:idx+20]!r}")
            node.name = parse_label()
            length = parse_length()
            return node, length

        label = parse_label()
        if not label:
            raise ValueError(f"Expected leaf label near {text[idx:idx+20]!r}")
        node = Node(label)
        length = parse_length()
        return node, length

    root, root_length = parse_subtree()
    if abs(root_length) > 0.0:
        raise ValueError("Root branch lengths are not supported")
    skip_ws()
    if idx != len(text):
        raise ValueError(f"Unexpected trailing Newick text: {text[idx:]!r}")
    return root


def all_nodes(root):
    seen = set()
    stack = [root]
    nodes = []
    while stack:
        node = stack.pop()
        if node.id in seen:
            continue
        seen.add(node.id)
        nodes.append(node)
        for other, _ in node.neighbors:
            stack.append(other)
    return nodes


def leaves_away_from(node, parent):
    leaves = []

    def visit(current, dad):
        if current.name:
            leaves.append(current.name)
        for child, _ in current.neighbors:
            if child is dad:
                continue
            visit(child, current)

    visit(node, parent)
    return set(leaves)


def find_edge_for_split(root, target_taxa):
    target = set(target_taxa)
    taxa = {node.name for node in all_nodes(root) if node.name}
    for node in all_nodes(root):
        for other, length in node.neighbors:
            if node.id > other.id:
                continue
            side = leaves_away_from(node, other)
            other_side = taxa - side
            if side == target or other_side == target:
                return node, other, length, side, other_side
    return None


def parse_fasta(path):
    sequences = {}
    name = None
    chunks = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        sequences[name] = "".join(chunks).upper()
    return sequences


def parse_phylip(path):
    with open(path, "r", encoding="utf-8") as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    if not lines:
        return {}
    header = lines[0].split()
    if len(header) < 2:
        raise ValueError(f"Invalid PHYLIP header in {path}")
    ntax = int(header[0])
    nsites = int(header[1])
    sequences = {}
    line_index = 1
    for _ in range(ntax):
        if line_index >= len(lines):
            raise ValueError(f"Unexpected end of PHYLIP file in {path}")
        parts = lines[line_index].split()
        line_index += 1
        if len(parts) < 2:
            raise ValueError(f"Invalid PHYLIP sequence row in {path}: {lines[line_index - 1]!r}")
        name = parts[0]
        chunks = ["".join(parts[1:])]
        while len("".join(chunks)) < nsites and line_index < len(lines):
            chunks.append("".join(lines[line_index].split()))
            line_index += 1
        sequence = "".join(chunks).upper()
        if len(sequence) != nsites:
            raise ValueError(f"PHYLIP sequence length mismatch for {name} in {path}: expected {nsites}, got {len(sequence)}")
        sequences[name] = sequence
    return sequences


def parse_alignment(path):
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                return parse_fasta(path)
            return parse_phylip(path)
    return {}


def parse_sat_stat(path, target_taxa):
    target = ",".join(sorted(target_taxa))
    with open(path, "r", encoding="utf-8") as handle:
        header = None
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("#") or not line:
                continue
            fields = line.split("\t")
            if fields[0] == "ID":
                header = fields
                continue
            if not header or not fields[0].isdigit():
                continue
            row = dict(zip(header, fields))
            if row.get("Formula", "dominant") != "dominant":
                continue
            if row.get("RateCategory", "pooled") != "pooled":
                continue
            split = row.get("Split", "")
            if split == target or set(split.split(",")) == set(target_taxa):
                return row
    return None


def parse_model(model):
    if model == "JC":
        return np.ones(6, dtype=float), np.repeat(0.25, 4)

    if not (model.startswith("GTR{") and "}+F{" in model and model.endswith("}")):
        raise ValueError(f"Unsupported model syntax: {model}")
    rate_text, freq_text = model[4:-1].split("}+F{", 1)
    rates = np.array([float(value) for value in rate_text.split(",")], dtype=float)
    pi = np.array([float(value) for value in freq_text.split(",")], dtype=float)
    if len(rates) != 6 or len(pi) != 4:
        raise ValueError(f"Unsupported model dimensions: {model}")
    pi /= pi.sum()
    return rates, pi


def build_q(model):
    rates, pi = parse_model(model)
    q = np.zeros((4, 4), dtype=float)
    for i, j, rate in [
        (0, 1, rates[0]),
        (0, 2, rates[1]),
        (0, 3, rates[2]),
        (1, 2, rates[3]),
        (1, 3, rates[4]),
        (2, 3, rates[5]),
    ]:
        q[i, j] = rate * pi[j]
        q[j, i] = rate * pi[i]
    for i in range(4):
        q[i, i] = -float(np.sum(q[i, :]))
    q /= -float(np.dot(pi, np.diag(q)))
    return q, pi


def reversible_eigendecomposition(q, pi):
    sqrt_pi = np.sqrt(pi)
    inv_sqrt_pi = 1.0 / sqrt_pi
    sym_q = np.diag(sqrt_pi) @ q @ np.diag(inv_sqrt_pi)
    evals, basis = np.linalg.eigh(sym_q)
    eigenvectors = np.diag(inv_sqrt_pi) @ basis
    inv_eigenvectors = basis.T @ np.diag(sqrt_pi)
    return evals, eigenvectors, inv_eigenvectors


def transition_matrix(evals, eigenvectors, inv_eigenvectors, length):
    exp_diag = np.diag(np.exp(evals * length))
    transition = eigenvectors @ exp_diag @ inv_eigenvectors
    transition[np.logical_and(transition < 0.0, transition > -1e-14)] = 0.0
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


def branch_key(a, b):
    return tuple(sorted((a.id, b.id)))


def compute_partial(node, parent, site, sequences, transitions, memo):
    key = (node.id, -1 if parent is None else parent.id, site)
    if key in memo:
        return memo[key]

    if node.name:
        symbol = sequences[node.name][site]
        state = STATE_INDEX.get(symbol)
        if state is None:
            partial = np.ones(4, dtype=float)
        else:
            partial = np.zeros(4, dtype=float)
            partial[state] = 1.0
        memo[key] = partial
        return partial

    partial = np.ones(4, dtype=float)
    for child, length in node.neighbors:
        if child is parent:
            continue
        child_partial = compute_partial(child, node, site, sequences, transitions, memo)
        contribution = transitions[branch_key(node, child)] @ child_partial
        partial *= contribution
    memo[key] = partial
    return partial


def formula_weights(evals, modes, formula, branch_length):
    if formula in {"dominant", "all_unweighted"}:
        return np.ones(len(modes), dtype=float)
    dominant = max(evals[m] for m in modes)
    return np.array([math.exp((evals[m] - dominant) * branch_length) for m in modes], dtype=float)


def compute_formula(root, sequences, target_taxa, model, formula, alpha, alpha_used):
    results = compute_formulas(root, sequences, target_taxa, model, [formula], alpha, alpha_used)
    return results.get(formula)


def compute_formulas(root, sequences, target_taxa, model, formulas, alpha, alpha_used):
    edge = find_edge_for_split(root, target_taxa)
    if edge is None:
        return {formula: None for formula in formulas}
    left, right, branch_length, left_taxa, right_taxa = edge
    q, pi = build_q(model)
    evals, eigenvectors, inv_eigenvectors = reversible_eigendecomposition(q, pi)
    transitions = {
        branch_key(node, other): transition_matrix(evals, eigenvectors, inv_eigenvectors, length)
        for node in all_nodes(root)
        for other, length in node.neighbors
        if node.id < other.id
    }
    formula_modes = {formula: mode_indices(evals, formula) for formula in formulas}
    formula_weights_map = {
        formula: formula_weights(evals, formula_modes[formula], formula, branch_length)
        for formula in formulas
    }
    formula_state = {}
    for formula, modes in formula_modes.items():
        formula_state[formula] = {
            "coherence": 0.0,
            "left_second": np.zeros((len(modes), len(modes)), dtype=float),
            "right_second": np.zeros((len(modes), len(modes)), dtype=float),
            "valid_sites": 0,
            "skipped_sites": 0,
        }

    nsites = len(next(iter(sequences.values())))

    for site in range(nsites):
        memo = {}
        left_lh = compute_partial(left, right, site, sequences, transitions, memo)
        right_lh = compute_partial(right, left, site, sequences, transitions, memo)
        left_post = left_lh * pi
        right_post = right_lh * pi
        left_sum = float(np.sum(left_post))
        right_sum = float(np.sum(right_post))
        if left_sum <= 0.0 or right_sum <= 0.0 or not math.isfinite(left_sum + right_sum):
            for state in formula_state.values():
                state["skipped_sites"] += 1
            continue
        left_post /= left_sum
        right_post /= right_sum

        full_left_factor = np.array([float(np.dot(eigenvectors[:, mode], left_post)) for mode in range(len(evals))], dtype=float)
        full_right_factor = np.array([float(np.dot(eigenvectors[:, mode], right_post)) for mode in range(len(evals))], dtype=float)
        if not np.all(np.isfinite(full_left_factor)) or not np.all(np.isfinite(full_right_factor)):
            for state in formula_state.values():
                state["skipped_sites"] += 1
            continue

        for formula, state in formula_state.items():
            modes = formula_modes[formula]
            weights = formula_weights_map[formula]
            left_factor = full_left_factor[modes]
            right_factor = full_right_factor[modes]
            state["coherence"] += float(np.dot(weights, left_factor * right_factor))
            state["left_second"] += np.outer(left_factor, left_factor)
            state["right_second"] += np.outer(right_factor, right_factor)
            state["valid_sites"] += 1

    results = {}
    for formula, state in formula_state.items():
        modes = formula_modes[formula]
        weights = formula_weights_map[formula]
        valid_sites = state["valid_sites"]
        if valid_sites == 0:
            results[formula] = None
            continue

        coherence = state["coherence"] / valid_sites
        weight_matrix = np.outer(weights, weights)
        variance = float(
            np.sum(weight_matrix * (state["left_second"] / valid_sites) * (state["right_second"] / valid_sites))
        )
        se = math.sqrt(variance / valid_sites) if variance > 0.0 else math.nan
        z_score = coherence / se if se > 0.0 and math.isfinite(se) else math.nan
        p_value = 0.5 * math.erfc(z_score / math.sqrt(2.0)) if math.isfinite(z_score) else math.nan
        results[formula] = {
            "target_found": 1,
            "left_taxa": len(left_taxa),
            "right_taxa": len(right_taxa),
            "valid_sites": valid_sites,
            "skipped_sites": state["skipped_sites"],
            "branch_length_used": branch_length,
            "satC": coherence,
            "satVar": variance,
            "satSE": se,
            "satZ": z_score,
            "satP": p_value,
            "decision": "informative" if math.isfinite(p_value) and p_value <= alpha_used else "saturated",
            "alpha": alpha,
            "alpha_used": alpha_used,
            "modes": ",".join(str(m) for m in modes),
            "eigenvalues": ",".join(f"{evals[m]:.10g}" for m in modes),
            "weights": ",".join(f"{w:.10g}" for w in weights),
        }
    return results


def pair(prefix, i, j, external_length, stem_length):
    return f"({prefix}{i}:{external_length:.6g},{prefix}{j}:{external_length:.6g}):{stem_length:.6g}"


def balanced_four_taxon_subtree(prefix, length=0.2):
    return f"(({prefix}1:{length:.6g},{prefix}2:{length:.6g}):{length:.6g},({prefix}3:{length:.6g},{prefix}4:{length:.6g}):{length:.6g})"


def eight_taxon_subtree(prefix, rng, pools):
    def ext():
        return sample_length(rng, pools, "external", 0.1)

    def internal():
        return sample_length(rng, pools, "internal", 0.1)

    left = f"({pair(prefix, 1, 2, ext(), internal())},{pair(prefix, 3, 4, ext(), internal())}):{internal():.6g}"
    right = f"({pair(prefix, 5, 6, ext(), internal())},{pair(prefix, 7, 8, ext(), internal())}):{internal():.6g}"
    return f"({left},{right})"


def five_taxon_tree(branch_length):
    subtree_a = balanced_four_taxon_subtree("A", 0.2)
    return f"({subtree_a}:0.0,B:{branch_length:.6g});\n"


def sixteen_taxon_tree(branch_length, rng, pools):
    return f"({eight_taxon_subtree('A', rng, pools)}:{branch_length:.6g},{eight_taxon_subtree('B', rng, pools)}:0.0);\n"


def target_taxa_for_case(tree_case):
    if tree_case == "five_external":
        return ["B"]
    if tree_case == "sixteen_internal":
        return [f"A{i}" for i in range(1, 9)]
    raise ValueError(tree_case)


def write_tree(tree_case, branch_length, path, rng, pools):
    if tree_case == "five_external":
        path.write_text(five_taxon_tree(branch_length), encoding="utf-8")
    elif tree_case == "sixteen_internal":
        path.write_text(sixteen_taxon_tree(branch_length, rng, pools), encoding="utf-8")
    else:
        raise ValueError(tree_case)


def run_seqgen(seqgen, tree_file, model, nsites, seed, alignment):
    if not seqgen or not os.path.exists(seqgen):
        raise SystemExit(
            "Seq-Gen is required for a paper-faithful rerun but was not found. "
            "Pass --seqgen /path/to/seq-gen, or use --simulator alisim for smoke tests."
        )
    if model == "JC":
        cmd = [seqgen, "-mHKY", "-l", str(nsites), "-n", "1", "-z", str(seed)]
    elif model.startswith("GTR{"):
        rates, pi = parse_model(model)
        cmd = [
            seqgen,
            "-mGTR",
            "-r",
            ",".join(f"{value:.10g}" for value in rates),
            "-f",
            ",".join(f"{value:.10g}" for value in pi),
            "-l",
            str(nsites),
            "-n",
            "1",
            "-z",
            str(seed),
        ]
    else:
        raise ValueError(f"Seq-Gen simulator is not configured for model {model}")

    with open(tree_file, "r", encoding="utf-8") as tree_handle, open(alignment, "w", encoding="utf-8") as out_handle:
        subprocess.run(cmd, stdin=tree_handle, stdout=out_handle, stderr=subprocess.DEVNULL, check=True)


def simulate_alignment(iqtree, seqgen, simulator, tree_file, model, nsites, seed, sim_prefix):
    if simulator == "seq-gen":
        alignment = Path(str(sim_prefix) + ".phy")
        run_seqgen(seqgen, tree_file, model, nsites, seed, alignment)
        return alignment

    run([iqtree, "--alisim", str(sim_prefix), "-t", str(tree_file), "-m", model, "--length", str(nsites), "--seed", str(seed), "-af", "fasta", "--quiet"])
    return Path(str(sim_prefix) + ".fa")


def run_satute(iqtree, alignment, prefix, model, tree=None, fixed_lengths=False):
    cmd = [iqtree, "-s", str(alignment)]
    if tree is not None:
        cmd.extend(["-te", str(tree)])
    cmd.extend(["-m", model, "--satute", "--prefix", str(prefix), "-T", "1", "--redo", "--quiet"])
    if fixed_lengths:
        cmd.append("-blfix")
    run(cmd)
    return Path(str(prefix) + ".sat.tree"), Path(str(prefix) + ".sat.stat")


def write_missing_rows(writer, base, scenario, formula, alpha_used):
    row = dict(base)
    row.update(
        {
            "scenario": scenario,
            "formula": formula,
            "target_found": 0,
            "left_taxa": "",
            "right_taxa": "",
            "valid_sites": "",
            "skipped_sites": "",
            "branch_length_used": "",
            "alpha_used": alpha_used,
            "satC": "",
            "satVar": "",
            "satSE": "",
            "satZ": "",
            "satP": "",
            "decision": "",
            "modes": "",
            "eigenvalues": "",
            "weights": "",
            "iqtree_satZ": "",
            "iqtree_satP": "",
            "iqtree_decision": "",
            "iqtree_bonf_decision": "",
        }
    )
    writer.writerow(row)


def write_formula_rows(writer, base, scenario, formulas, tree_path, stat_path, alignment, model, target_taxa, alpha, alpha_used):
    sequences = parse_alignment(alignment)
    root = parse_newick(tree_path)
    iqtree_row = parse_sat_stat(stat_path, target_taxa)
    results = compute_formulas(root, sequences, target_taxa, model, formulas, alpha, alpha_used)
    for formula in formulas:
        result = results.get(formula)
        if result is None:
            write_missing_rows(writer, base, scenario, formula, alpha_used)
            continue
        row = dict(base)
        row.update({"scenario": scenario, "formula": formula})
        row.update(result)
        row["iqtree_satZ"] = iqtree_row.get("satZ", "") if iqtree_row else ""
        row["iqtree_satP"] = iqtree_row.get("satP", "") if iqtree_row else ""
        row["iqtree_decision"] = iqtree_row.get("Decision", "") if iqtree_row else ""
        row["iqtree_bonf_decision"] = iqtree_row.get("DecisionTaxonBonf", "") if iqtree_row else ""
        writer.writerow(row)


def run_case(iqtree, seqgen, simulator, pools, outdir, writer, tree_case, model, nsites, branch_length, rep, seed, alpha):
    case_dir = outdir / "runs" / tree_case / sanitize_model(model) / f"n{nsites}" / f"b{branch_length:.2f}" / f"r{rep:04d}"
    case_dir.mkdir(parents=True, exist_ok=True)
    tree_file = case_dir / "true.tree"
    rng = random.Random(seed)
    write_tree(tree_case, branch_length, tree_file, rng, pools)

    sim_prefix = case_dir / "sim"
    alignment = simulate_alignment(iqtree, seqgen, simulator, tree_file, model, nsites, seed, sim_prefix)
    target_taxa = target_taxa_for_case(tree_case)
    formulas = ["dominant", "all_unweighted", "eigenvalue_weighted"]
    base = {
        "tree_case": tree_case,
        "model": model,
        "nsites": nsites,
        "branch_length": branch_length,
        "replicate": rep,
        "seed": seed,
        "simulator": simulator,
        "branch_length_source": "evonaps_table" if pools.get("internal") and pools.get("external") else "fixed_0.1_smoke",
        "target_split": ",".join(target_taxa),
        "alpha": alpha,
    }

    tree_path, stat_path = run_satute(iqtree, alignment, case_dir / "true_fixed", model, tree_file, True)
    write_formula_rows(writer, base, "true_tree_fixed_lengths", formulas, tree_path, stat_path, alignment, model, target_taxa, alpha, alpha)

    tree_path, stat_path = run_satute(iqtree, alignment, case_dir / "true_ml_lengths", model, tree_file, False)
    write_formula_rows(writer, base, "true_topology_ml_lengths", formulas, tree_path, stat_path, alignment, model, target_taxa, alpha, alpha)

    tree_path, stat_path = run_satute(iqtree, alignment, case_dir / "ml_tree", model, None, False)
    write_formula_rows(writer, base, "ml_tree_unadjusted", formulas, tree_path, stat_path, alignment, model, target_taxa, alpha, alpha)

    taxa_count = 5 if tree_case == "five_external" else 16
    alpha_bonf = alpha / (1 * (taxa_count - 1) if tree_case == "five_external" else 8 * 8)
    write_formula_rows(writer, base, "ml_tree_bonferroni", formulas, tree_path, stat_path, alignment, model, target_taxa, alpha, alpha_bonf)


def aggregate(detail_path, summary_path):
    groups = defaultdict(lambda: {"evaluated": 0, "informative": 0, "missing": 0})
    with open(detail_path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = (row["tree_case"], row["model"], row["nsites"], row["branch_length"], row["scenario"], row["formula"])
            if row["target_found"] == "1":
                groups[key]["evaluated"] += 1
                groups[key]["informative"] += 1 if row["decision"] == "informative" else 0
            else:
                groups[key]["missing"] += 1

    with open(summary_path, "w", encoding="utf-8", newline="") as handle:
        fieldnames = ["tree_case", "model", "nsites", "branch_length", "scenario", "formula", "evaluated", "informative", "missing_split", "fraction_informative"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for key in sorted(groups, key=lambda x: (x[0], x[1], int(x[2]), float(x[3]), x[4], x[5])):
            group = groups[key]
            fraction = group["informative"] / group["evaluated"] if group["evaluated"] else ""
            writer.writerow(
                {
                    "tree_case": key[0],
                    "model": key[1],
                    "nsites": key[2],
                    "branch_length": key[3],
                    "scenario": key[4],
                    "formula": key[5],
                    "evaluated": group["evaluated"],
                    "informative": group["informative"],
                    "missing_split": group["missing"],
                    "fraction_informative": fraction,
                }
            )


def sanitize_model(model):
    return model.replace("{", "_").replace("}", "_").replace("+", "_").replace(",", "_").replace(".", "p")


def task_key_from_row(row):
    return (
        row["tree_case"],
        row["model"],
        str(row["nsites"]),
        f"{float(row['branch_length']):.10g}",
        str(row["replicate"]),
    )


def completed_tasks(detail_path, fieldnames=None, clean_incomplete=False):
    if not detail_path.exists():
        return set()
    rows = []
    counts = defaultdict(int)
    with open(detail_path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append(row)
            key = task_key_from_row(row)
            counts[key] += 1
    complete = {key for key, count in counts.items() if count >= 12}
    if clean_incomplete and fieldnames is not None:
        with open(detail_path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for row in rows:
                if task_key_from_row(row) in complete:
                    writer.writerow(row)
    return complete


def main():
    parser = argparse.ArgumentParser(description="Paired head-to-head SatuTe formula simulations.")
    parser.add_argument("--iqtree", default="./build/iqtree3")
    parser.add_argument("--outdir", default="/tmp/iqtree-satute-head-to-head")
    parser.add_argument("--reps", type=int, default=1)
    parser.add_argument("--site-lengths", default="100")
    parser.add_argument("--branch-lengths", default="0.1,1.0")
    parser.add_argument("--tree-cases", default="five_external")
    parser.add_argument("--models", default="JC", help="Comma-separated model aliases, e.g. JC or GTR_PF06346.")
    parser.add_argument("--paper-grid", action="store_true", help="Use paper branch lengths and site lengths 100,1000,10000.")
    parser.add_argument("--simulator", choices=["alisim", "seq-gen"], default="alisim")
    parser.add_argument("--seqgen", default=shutil.which("seq-gen") or shutil.which("seqgen") or "")
    parser.add_argument("--evonaps-branch-lengths", default="", help="CSV/TSV with columns type/kind/branch_type and length/branch_length/blen.")
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--resume", action="store_true", help="Append to an existing shard output and skip completed replicate tasks.")
    args = parser.parse_args()

    iqtree = str(Path(args.iqtree))
    if not os.path.exists(iqtree) or not os.access(iqtree, os.X_OK):
        raise SystemExit(f"Cannot execute IQ-TREE binary: {iqtree}")

    outdir = Path(args.outdir)
    if outdir.exists() and not args.resume:
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if args.shard_count < 1:
        raise SystemExit("--shard-count must be >= 1")
    if args.shard_index < 0 or args.shard_index >= args.shard_count:
        raise SystemExit("--shard-index must satisfy 0 <= shard-index < shard-count")

    branch_lengths = parse_csv_numbers(PAPER_BRANCH_LENGTHS if args.paper_grid else args.branch_lengths, float)
    site_lengths = parse_csv_numbers("100,1000,10000" if args.paper_grid else args.site_lengths, int)
    tree_cases = [value.strip() for value in args.tree_cases.split(",") if value.strip()]
    model_aliases = [value.strip() for value in args.models.split(",") if value.strip()]
    models = [resolve_model_alias(alias) for alias in model_aliases]
    pools = load_branch_length_pools(args.evonaps_branch_lengths)

    detail_path = outdir / "head_to_head_detail.tsv"
    summary_path = outdir / "head_to_head_summary.tsv"
    fieldnames = [
        "tree_case",
        "model",
        "nsites",
        "branch_length",
        "replicate",
        "seed",
        "simulator",
        "branch_length_source",
        "target_split",
        "scenario",
        "formula",
        "target_found",
        "left_taxa",
        "right_taxa",
        "valid_sites",
        "skipped_sites",
        "branch_length_used",
        "alpha",
        "alpha_used",
        "satC",
        "satVar",
        "satSE",
        "satZ",
        "satP",
        "decision",
        "modes",
        "eigenvalues",
        "weights",
        "iqtree_satZ",
        "iqtree_satP",
        "iqtree_decision",
        "iqtree_bonf_decision",
    ]
    done = completed_tasks(detail_path, fieldnames, clean_incomplete=True) if args.resume else set()
    write_header = not (args.resume and detail_path.exists() and detail_path.stat().st_size > 0)
    mode = "a" if args.resume else "w"
    with open(detail_path, mode, encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        if write_header:
            writer.writeheader()
        task_index = 0
        selected_tasks = 0
        skipped_tasks = 0
        for tree_case in tree_cases:
            for model in models:
                for nsites in site_lengths:
                    for branch_length in branch_lengths:
                        for rep in range(1, args.reps + 1):
                            current = task_index
                            task_index += 1
                            if current % args.shard_count != args.shard_index:
                                continue
                            selected_tasks += 1
                            task_key = (tree_case, model, str(nsites), f"{branch_length:.10g}", str(rep))
                            if task_key in done:
                                skipped_tasks += 1
                                continue
                            seed = 900000 + rep + nsites * 10 + int(branch_length * 1000)
                            run_case(iqtree, args.seqgen, args.simulator, pools, outdir, writer, tree_case, model, nsites, branch_length, rep, seed, args.alpha)

    aggregate(detail_path, summary_path)
    print(f"Shard:   {args.shard_index}/{args.shard_count}")
    print(f"Tasks:   {selected_tasks}/{task_index}")
    if args.resume:
        print(f"Skipped: {skipped_tasks}")
    print(f"Detail:  {detail_path}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
