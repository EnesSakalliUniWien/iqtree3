#!/usr/bin/env python3

import argparse
import csv
from contextlib import ExitStack
import hashlib
import math
import os
import random
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
IQTREE_DEFAULT = PROJECT_ROOT.parents[1] / "build" / "iqtree3"
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from satute_analysis.benchmark_contracts import (  # noqa: E402
    BASE_DETAIL_FIELDS,
    BRANCH_AUDIT_FIELDS,
    FDR_AUDIT_FIELDS,
    FORMULAS,
    ML_TREE,
    SCHEMA_VERSION,
    SUMMARY_FIELDS,
    TIMING_FIELDS,
    TRUE_FIXED,
    TRUE_ML_LENGTHS,
    expected_rows_for_task,
    selected_scenarios,
)
from satute_analysis.benchmark_summary import (  # noqa: E402
    aggregate_detail,
    validate_tsv_header,
)
from satute_analysis.native_results import (  # noqa: E402
    branch_audit_detail,
    fdr_family_audit,
    index_target_rows,
    read_native_rows,
    target_detail,
)

# The production benchmark consumes native IQ-TREE rows and never imports
# NumPy. The dormant independent-reference helper loads it only when called.
np = None


DNA_STATE_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3, "U": 3}
AA_ALPHABET = "ARNDCQEGHILKMFPSTWYV"
AA_STATE_INDEX = {symbol: index for index, symbol in enumerate(AA_ALPHABET)}
PAPER_BRANCH_LENGTHS = "0.1,0.2,0.3,0.4,0.5,0.8,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,7.5,10.0"
GTR_PF06346_MODEL = "GTR{0.6676,3.7807,4.2833,0.5354,0.8718,1.0}+F{0.125,0.436,0.191,0.245}"
GTR_SKEW_FREQ_MODEL = "GTR{1.0,1.0,1.0,1.0,1.0,1.0}+F{0.70,0.10,0.10,0.10}"
GTR_SKEW_RATES_MODEL = "GTR{0.05,8.0,0.10,0.10,5.0,0.05}+F{0.25,0.25,0.25,0.25}"
GTR_SKEW_BOTH_MODEL = "GTR{0.05,8.0,0.10,0.10,5.0,0.05}+F{0.70,0.10,0.10,0.10}"
PROTEIN_GAMMA_SHAPE = 0.5
LG_RATES = [
    0.425093, 0.276818, 0.395144, 2.489084, 0.969894, 1.038545, 2.066040, 0.358858, 0.149830,
    0.395337, 0.536518, 1.124035, 0.253701, 1.177651, 4.727182, 2.139501, 0.180717, 0.218959,
    2.547870, 0.751878, 0.123954, 0.534551, 2.807908, 0.363970, 0.390192, 2.426601, 0.126991,
    0.301848, 6.326067, 0.484133, 0.052722, 0.332533, 0.858151, 0.578987, 0.593607, 0.314440,
    0.170887, 5.076149, 0.528768, 1.695752, 0.541712, 1.437645, 4.509238, 0.191503, 0.068427,
    2.145078, 0.371004, 0.089525, 0.161787, 4.008358, 2.000679, 0.045376, 0.612025, 0.083688,
    0.062556, 0.523386, 5.243870, 0.844926, 0.927114, 0.010690, 0.015076, 0.282959, 0.025548,
    0.017416, 0.394456, 1.240275, 0.425860, 0.029890, 0.135107, 0.037967, 0.084808, 0.003499,
    0.569265, 0.640543, 0.320627, 0.594007, 0.013266, 0.893680, 1.105251, 0.075382, 2.784478,
    1.143480, 0.670128, 1.165532, 1.959291, 4.128591, 0.267959, 4.813505, 0.072854, 0.582457,
    3.234294, 1.672569, 0.035855, 0.624294, 1.223828, 1.080136, 0.236199, 0.257336, 0.210332,
    0.348847, 0.423881, 0.044265, 0.069673, 1.807177, 0.173735, 0.018811, 0.419409, 0.611973,
    0.604545, 0.077852, 0.120037, 0.245034, 0.311484, 0.008705, 0.044261, 0.296636, 0.139538,
    0.089586, 0.196961, 1.739990, 0.129836, 0.268491, 0.054679, 0.076701, 0.108882, 0.366317,
    0.697264, 0.442472, 0.682139, 0.508851, 0.990012, 0.584262, 0.597054, 5.306834, 0.119013,
    4.145067, 0.159069, 4.273607, 1.112727, 0.078281, 0.064105, 1.033739, 0.111660, 0.232523,
    10.649107, 0.137500, 6.312358, 2.592692, 0.249060, 0.182287, 0.302936, 0.619632, 0.299648,
    1.702745, 0.656604, 0.023918, 0.390322, 0.748683, 1.136863, 0.049906, 0.131932, 0.185202,
    1.798853, 0.099849, 0.346960, 2.020366, 0.696175, 0.481306, 1.898718, 0.094464, 0.361819,
    0.165001, 2.457121, 7.803902, 0.654683, 1.338132, 0.571468, 0.095131, 0.089613, 0.296501,
    6.472279, 0.248862, 0.400547, 0.098369, 0.140825, 0.245841, 2.188158, 3.151815, 0.189510,
    0.249313,
]
LG_FREQ = [
    0.079066, 0.055941, 0.041977, 0.053052, 0.012937, 0.040767, 0.071586, 0.057337, 0.022355,
    0.062157, 0.099081, 0.064600, 0.022951, 0.042302, 0.044040, 0.061197, 0.053287, 0.012066,
    0.034155, 0.069147,
]
MODEL_ALIASES = {
    "JC": "JC",
    "LG": "LG",
    "WAG": "WAG",
    "JTT": "JTT",
    "Q.PFAM": "Q.pfam",
    "Q.pfam": "Q.pfam",
    "LG_G4": f"LG+G4{{{PROTEIN_GAMMA_SHAPE}}}",
    "WAG_G4": f"WAG+G4{{{PROTEIN_GAMMA_SHAPE}}}",
    "JTT_G4": f"JTT+G4{{{PROTEIN_GAMMA_SHAPE}}}",
    "Q.PFAM_G4": f"Q.pfam+G4{{{PROTEIN_GAMMA_SHAPE}}}",
    "Q.pfam_G4": f"Q.pfam+G4{{{PROTEIN_GAMMA_SHAPE}}}",
    "K2P": "K2P",
    "F81": "F81",
    "GTR_PF06346": GTR_PF06346_MODEL,
    "GTR_EvoNAPS_PF06346": GTR_PF06346_MODEL,
    "GTR_SKEW_FREQ": GTR_SKEW_FREQ_MODEL,
    "GTR_SKEW_RATES": GTR_SKEW_RATES_MODEL,
    "GTR_SKEW_BOTH": GTR_SKEW_BOTH_MODEL,
}
FIXED_EMPIRICAL_PROTEIN_MODELS = {"LG", "WAG", "JTT", "Q.pfam"}
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


def parse_csv_text(text):
    return [value.strip() for value in text.split(",") if value.strip()]


def simulation_seed(base_seed, tree_case, simulation_model, nsites, branch_length, replicate):
    """Return a deterministic positive seed shared across evaluation models."""

    payload = (
        f"{base_seed}:{tree_case}:{simulation_model}:{nsites}:"
        f"{float(branch_length):.12g}:{replicate}"
    ).encode("utf-8")
    value = int.from_bytes(hashlib.blake2s(payload, digest_size=8).digest(), "big")
    return 1 + value % 2147483646


def resolve_model_alias(alias):
    if alias in MODEL_ALIASES:
        return MODEL_ALIASES[alias]
    return alias


def split_fixed_gamma(model):
    """Return the substitution model and an explicitly fixed G4 shape."""
    marker = "+G4{"
    if marker not in model:
        return model, None
    base, suffix = model.rsplit(marker, 1)
    if not suffix.endswith("}"):
        raise ValueError(f"Invalid fixed gamma model syntax: {model}")
    alpha = float(suffix[:-1])
    if not math.isfinite(alpha) or alpha <= 0.0:
        raise ValueError(f"Gamma shape must be positive and finite: {model}")
    return base, alpha


def is_fixed_model_spec(model):
    base, _gamma_shape = split_fixed_gamma(resolve_model_alias(model))
    return (
        base in {"JC", *FIXED_EMPIRICAL_PROTEIN_MODELS}
        or base.startswith("GTR{")
        or base.startswith("GTR20{")
        or base.startswith("F81+F{")
        or base.startswith("K2P{")
    )


def parse_model_pairs(text):
    pairs = []
    for item in parse_csv_text(text):
        if ":" not in item:
            raise ValueError(f"Model pair must have SIM:EVAL syntax, got: {item}")
        simulation_alias, evaluation_alias = [part.strip() for part in item.split(":", 1)]
        if not simulation_alias or not evaluation_alias:
            raise ValueError(f"Model pair must have non-empty SIM:EVAL aliases, got: {item}")
        pairs.append((resolve_model_alias(simulation_alias), resolve_model_alias(evaluation_alias)))
    return pairs


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
    return index_target_rows(read_native_rows(path), target_taxa)


def parse_model(model):
    model, _gamma_shape = split_fixed_gamma(model)
    if model == "LG":
        return list(LG_RATES), list(LG_FREQ)

    if model == "JC":
        return [1.0] * 6, [0.25] * 4

    if model.startswith("F81+F{") and model.endswith("}"):
        freq_text = model[len("F81+F{") : -1]
        pi = [float(value) for value in freq_text.split(",")]
        if len(pi) != 4:
            raise ValueError(f"Unsupported F81 model dimensions: {model}")
        total = sum(pi)
        pi = [value / total for value in pi]
        return [1.0] * 6, pi

    if model.startswith("K2P{") and model.endswith("}+FQ"):
        kappa = float(model[len("K2P{") : -len("}+FQ")])
        return [1.0, kappa, 1.0, 1.0, kappa, 1.0], [0.25] * 4

    if model in {"F81", "K2P"}:
        raise ValueError(f"Model {model} requires fitted parameters from IQ-TREE before independent calculation")

    is_gtr4 = model.startswith("GTR{")
    is_gtr20 = model.startswith("GTR20{")
    if not ((is_gtr4 or is_gtr20) and "}+F{" in model and model.endswith("}")):
        raise ValueError(f"Unsupported model syntax: {model}")
    prefix_length = len("GTR20{") if is_gtr20 else len("GTR{")
    rate_text, freq_text = model[prefix_length:-1].split("}+F{", 1)
    rates = [float(value) for value in rate_text.split(",")]
    pi = [float(value) for value in freq_text.split(",")]
    expected_rate_count = len(pi) * (len(pi) - 1) // 2
    # IQ-TREE fixes the final GTR20 exchangeability to one for scale and accepts
    # 189 values. Seq-Gen GENERAL requires the complete 190-value triangle.
    if len(pi) == 20 and len(rates) == expected_rate_count - 1:
        rates.append(1.0)
    if len(pi) not in {4, 20} or len(rates) != expected_rate_count:
        raise ValueError(f"Unsupported model dimensions: {model}")
    total = sum(pi)
    pi = [value / total for value in pi]
    return rates, pi


def build_q(model):
    rates, pi = parse_model(model)
    rates = np.asarray(rates, dtype=float)
    pi = np.asarray(pi, dtype=float)
    nstates = len(pi)
    expected_rate_count = nstates * (nstates - 1) // 2
    if len(rates) != expected_rate_count:
        raise ValueError(f"Rate count {len(rates)} does not match {nstates} states for model {model}")
    q = np.zeros((nstates, nstates), dtype=float)
    rate_index = 0
    for i in range(nstates - 1):
        for j in range(i + 1, nstates):
            rate = rates[rate_index]
            rate_index += 1
            q[i, j] = rate * pi[j]
            q[j, i] = rate * pi[i]
    for i in range(nstates):
        q[i, i] = -float(np.sum(q[i, :]))
    q /= -float(np.dot(pi, np.diag(q)))
    return q, pi


def state_index_for_model(model):
    _, pi = parse_model(model)
    if len(pi) == 4:
        return DNA_STATE_INDEX, 4
    if len(pi) == 20:
        return AA_STATE_INDEX, 20
    raise ValueError(f"Unsupported state count for model {model}: {len(pi)}")


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
    if formula in {"eigenvalue_weighted", "eigenvalue_weighted_gls"}:
        return nonzero
    if formula != "dominant":
        raise ValueError(f"Unknown SatuTe formula: {formula}")

    dominant = max(evals[i] for i in nonzero)
    tol = max(1e-8, abs(dominant) * 1e-6)
    return [i for i in nonzero if abs(evals[i] - dominant) <= tol]


def branch_key(a, b):
    return tuple(sorted((a.id, b.id)))


def compute_partial(node, parent, site, sequences, transitions, memo, state_index, nstates):
    key = (node.id, -1 if parent is None else parent.id, site)
    if key in memo:
        return memo[key]

    if node.name:
        symbol = sequences[node.name][site]
        state = state_index.get(symbol)
        if state is None:
            partial = np.ones(nstates, dtype=float)
        else:
            partial = np.zeros(nstates, dtype=float)
            partial[state] = 1.0
        memo[key] = partial
        return partial

    partial = np.ones(nstates, dtype=float)
    for child, length in node.neighbors:
        if child is parent:
            continue
        child_partial = compute_partial(child, node, site, sequences, transitions, memo, state_index, nstates)
        contribution = transitions[branch_key(node, child)] @ child_partial
        partial *= contribution
    memo[key] = partial
    return partial


def formula_weights(evals, modes, formula, branch_length):
    if formula == "dominant":
        return np.ones(len(modes), dtype=float)
    dominant = max(evals[m] for m in modes)
    return np.array([math.exp((evals[m] - dominant) * branch_length) for m in modes], dtype=float)


def gls_weights(sigma, target):
    scale = float(np.trace(sigma)) / len(target) if len(target) else 1.0
    ridge = max(scale, 1.0) * 1e-10
    try:
        if np.linalg.cond(sigma) < 1e10:
            weights = np.linalg.solve(sigma, target)
        else:
            weights = np.linalg.solve(sigma + ridge * np.eye(len(target)), target)
    except np.linalg.LinAlgError:
        weights = np.linalg.pinv(sigma) @ target
    if not np.all(np.isfinite(weights)):
        weights = np.linalg.pinv(sigma + ridge * np.eye(len(target))) @ target
    dominant = weights[int(np.argmax(target))]
    if abs(dominant) > 1e-12 and math.isfinite(float(dominant)):
        weights = weights / dominant
    return weights


def empty_projection_stats():
    return {
        "projection_angle_sites": 0,
        "projection_cosine": math.nan,
        "projection_sine": math.nan,
        "projection_left_length": math.nan,
        "projection_right_length": math.nan,
        "projection_length_delta": math.nan,
        "projection_abs_length_delta": math.nan,
    }


def compute_formula(root, sequences, target_taxa, model, formula, alpha, alpha_used):
    results = compute_formulas(root, sequences, target_taxa, model, [formula], alpha, alpha_used)
    return results.get(formula)


def alignment_patterns(sequences):
    sequence_names = tuple(sequences)
    nsites = len(next(iter(sequences.values())))
    patterns = {}
    for site in range(nsites):
        pattern = tuple(sequences[name][site] for name in sequence_names)
        if pattern in patterns:
            representative, frequency = patterns[pattern]
            patterns[pattern] = (representative, frequency + 1)
        else:
            patterns[pattern] = (site, 1)
    return patterns.values()


def compute_formulas(root, sequences, target_taxa, model, formulas, alpha, alpha_used):
    global np
    if np is None:
        import importlib

        np = importlib.import_module("numpy")

    edge = find_edge_for_split(root, target_taxa)
    if edge is None:
        return {formula: None for formula in formulas}
    left, right, branch_length, left_taxa, right_taxa = edge
    q, pi = build_q(model)
    state_index, nstates = state_index_for_model(model)
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
            "mode_coherence": np.zeros(len(modes), dtype=float),
            "left_second": np.zeros((len(modes), len(modes)), dtype=float),
            "right_second": np.zeros((len(modes), len(modes)), dtype=float),
            "valid_sites": 0,
            "skipped_sites": 0,
            "angle_sites": 0,
            "cosine_sum": 0.0,
            "sine_sum": 0.0,
            "left_length_sum": 0.0,
            "right_length_sum": 0.0,
            "length_delta_sum": 0.0,
            "abs_length_delta_sum": 0.0,
        }

    for site, frequency in alignment_patterns(sequences):
        memo = {}
        left_lh = compute_partial(left, right, site, sequences, transitions, memo, state_index, nstates)
        right_lh = compute_partial(right, left, site, sequences, transitions, memo, state_index, nstates)
        left_post = left_lh * pi
        right_post = right_lh * pi
        left_sum = float(np.sum(left_post))
        right_sum = float(np.sum(right_post))
        if left_sum <= 0.0 or right_sum <= 0.0 or not math.isfinite(left_sum + right_sum):
            for state in formula_state.values():
                state["skipped_sites"] += frequency
            continue
        left_post /= left_sum
        right_post /= right_sum

        full_left_factor = np.array([float(np.dot(eigenvectors[:, mode], left_post)) for mode in range(len(evals))], dtype=float)
        full_right_factor = np.array([float(np.dot(eigenvectors[:, mode], right_post)) for mode in range(len(evals))], dtype=float)
        if not np.all(np.isfinite(full_left_factor)) or not np.all(np.isfinite(full_right_factor)):
            for state in formula_state.values():
                state["skipped_sites"] += frequency
            continue

        for formula, state in formula_state.items():
            modes = formula_modes[formula]
            weights = formula_weights_map[formula]
            left_factor = full_left_factor[modes]
            right_factor = full_right_factor[modes]
            weighted_left = np.sqrt(weights) * left_factor
            weighted_right = np.sqrt(weights) * right_factor
            site_coherence = float(np.dot(weighted_left, weighted_right))
            state["coherence"] += frequency * site_coherence
            state["mode_coherence"] += frequency * left_factor * right_factor
            state["left_second"] += frequency * np.outer(left_factor, left_factor)
            state["right_second"] += frequency * np.outer(right_factor, right_factor)
            state["valid_sites"] += frequency
            left_length = float(np.linalg.norm(weighted_left))
            right_length = float(np.linalg.norm(weighted_right))
            length_delta = left_length - right_length
            state["left_length_sum"] += frequency * left_length
            state["right_length_sum"] += frequency * right_length
            state["length_delta_sum"] += frequency * length_delta
            state["abs_length_delta_sum"] += frequency * abs(length_delta)
            if left_length > 0.0 and right_length > 0.0:
                cosine = max(-1.0, min(1.0, site_coherence / (left_length * right_length)))
                state["cosine_sum"] += frequency * cosine
                state["sine_sum"] += frequency * math.sqrt(max(0.0, 1.0 - cosine * cosine))
                state["angle_sites"] += frequency

    results = {}
    for formula, state in formula_state.items():
        modes = formula_modes[formula]
        weights = formula_weights_map[formula]
        valid_sites = state["valid_sites"]
        if valid_sites == 0:
            results[formula] = None
            continue

        mode_coherence = state["mode_coherence"] / valid_sites
        sigma = (state["left_second"] / valid_sites) * (state["right_second"] / valid_sites)
        if formula == "eigenvalue_weighted_gls":
            weights = gls_weights(sigma, weights)
        coherence = float(np.dot(weights, mode_coherence))
        variance = float(weights @ sigma @ weights)
        se = math.sqrt(variance / valid_sites) if variance > 0.0 else math.nan
        z_score = coherence / se if se > 0.0 and math.isfinite(se) else math.nan
        p_value = 0.5 * math.erfc(z_score / math.sqrt(2.0)) if math.isfinite(z_score) else math.nan
        projection_stats = empty_projection_stats()
        projection_stats.update(
            {
                "projection_left_length": state["left_length_sum"] / valid_sites,
                "projection_right_length": state["right_length_sum"] / valid_sites,
                "projection_length_delta": state["length_delta_sum"] / valid_sites,
                "projection_abs_length_delta": state["abs_length_delta_sum"] / valid_sites,
                "projection_angle_sites": state["angle_sites"],
            }
        )
        if state["angle_sites"]:
            projection_stats["projection_cosine"] = state["cosine_sum"] / state["angle_sites"]
            projection_stats["projection_sine"] = state["sine_sum"] / state["angle_sites"]
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
            **projection_stats,
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
    return f"({subtree_a[1:-1]},B:{branch_length:.6g});\n"


def sixteen_taxon_tree(branch_length, rng, pools):
    subtree_a = eight_taxon_subtree("A", rng, pools)
    subtree_b = eight_taxon_subtree("B", rng, pools)
    return f"({subtree_a[1:-1]},{subtree_b}:{branch_length:.6g});\n"


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


def seqgen_command(seqgen, model, nsites, seed):
    base_model, gamma_shape = split_fixed_gamma(model)
    gamma_args = ["-a", f"{gamma_shape:.10g}", "-g", "4"] if gamma_shape is not None else []
    if base_model == "JC":
        cmd = [seqgen, "-mHKY"]
    elif base_model in {"LG", "WAG", "JTT"}:
        cmd = [seqgen, f"-m{base_model}"]
    elif base_model.startswith("GTR{") or base_model.startswith("GTR20{"):
        rates, pi = parse_model(base_model)
        seqgen_model = "GTR" if len(pi) == 4 else "GENERAL"
        cmd = [
            seqgen,
            f"-m{seqgen_model}",
            "-r",
            ",".join(f"{value:.10g}" for value in rates),
            "-f",
            ",".join(f"{value:.10g}" for value in pi),
        ]
    else:
        raise ValueError(
            f"Seq-Gen simulator is not configured for model {model}; "
            "use --simulator alisim for IQ-TREE protein models such as Q.pfam."
        )
    return [*cmd, *gamma_args, "-l", str(nsites), "-n", "1", "-z", str(seed)]


def run_seqgen(seqgen, tree_file, model, nsites, seed, alignment):
    if not seqgen or not os.path.exists(seqgen):
        raise SystemExit(
            "Seq-Gen is required for a paper-faithful rerun but was not found. "
            "Pass --seqgen /path/to/seq-gen, or use --simulator alisim for smoke tests."
        )
    cmd = seqgen_command(seqgen, model, nsites, seed)

    with open(tree_file, "r", encoding="utf-8") as tree_handle, open(alignment, "w", encoding="utf-8") as out_handle:
        subprocess.run(cmd, stdin=tree_handle, stdout=out_handle, stderr=subprocess.DEVNULL, check=True)


def indelible_model_lines(model):
    rates, pi = parse_model(model)
    if len(pi) != 4:
        raise ValueError(f"INDELible backend currently supports nucleotide models only, got: {model}")
    if model == "JC":
        return ["[MODEL] m1", "  [submodel] JC"]
    # IQ-TREE writes DNA frequencies as A,C,G,T and GTR rates as
    # AC,AG,AT,CG,CT,GT. INDELible expects T,C,A,G frequencies and
    # CT,AT,GT,AC,CG,AG GTR exchangeabilities.
    indelible_rates = [rates[4], rates[2], rates[5], rates[0], rates[3], rates[1]]
    indelible_pi = [pi[3], pi[1], pi[0], pi[2]]
    return [
        "[MODEL] m1",
        "  [submodel] GTR " + " ".join(f"{value:.10g}" for value in indelible_rates),
        "  [statefreq] " + " ".join(f"{value:.10g}" for value in indelible_pi),
    ]


def run_indelible(indelible, tree_file, model, nsites, seed, sim_prefix):
    if not indelible or not os.path.exists(indelible):
        raise SystemExit(
            "INDELible is required for --simulator indelible but was not found. "
            "Pass --indelible /path/to/indelible."
        )
    workdir = Path(sim_prefix).parent
    output_name = Path(sim_prefix).name
    tree_text = Path(tree_file).read_text(encoding="utf-8").strip()
    if not tree_text.endswith(";"):
        tree_text += ";"
    control_lines = [
        "[TYPE] NUCLEOTIDE 1",
        "",
        "[SETTINGS]",
        "  [output] FASTA",
        f"  [randomseed] {seed}",
        "",
        *indelible_model_lines(model),
        "",
        f"[TREE] t1 {tree_text}",
        "",
        "[PARTITIONS] p1",
        f"  [t1 m1 {nsites}]",
        "",
        f"[EVOLVE] p1 1 {output_name}",
        "",
    ]
    (workdir / "control.txt").write_text("\n".join(control_lines), encoding="utf-8")
    subprocess.run([indelible], cwd=workdir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    alignment = workdir / f"{output_name}.fas"
    if not alignment.exists():
        raise FileNotFoundError(f"INDELible did not create expected FASTA output: {alignment}")
    return alignment


def simulate_alignment(iqtree, seqgen, indelible, simulator, tree_file, model, nsites, seed, sim_prefix):
    if simulator == "seq-gen":
        alignment = Path(str(sim_prefix) + ".phy")
        run_seqgen(seqgen, tree_file, model, nsites, seed, alignment)
        return alignment

    if simulator == "indelible":
        return run_indelible(indelible, tree_file, model, nsites, seed, sim_prefix)

    run([iqtree, "--alisim", str(sim_prefix), "-t", str(tree_file), "-m", model, "--length", str(nsites), "--seed", str(seed), "-af", "fasta", "--quiet"])
    return Path(str(sim_prefix) + ".fa")


def run_satute(iqtree, alignment, prefix, model, tree=None, fixed_lengths=False, seed=None):
    cmd = [iqtree, "-s", str(alignment)]
    if tree is not None:
        cmd.extend(["-te", str(tree)])
    cmd.extend(["-m", model, "--satute", "--prefix", str(prefix), "-T", "1", "--redo", "--quiet"])
    if seed is not None:
        cmd.extend(["--seed", str(seed)])
    if fixed_lengths:
        cmd.append("-blfix")
    run(cmd)
    return Path(str(prefix) + ".sat.tree"), Path(str(prefix) + ".sat.stat")


def fitted_model_for_compute(prefix, requested_model):
    requested_model = resolve_model_alias(requested_model)
    if is_fixed_model_spec(requested_model):
        return requested_model

    iqtree_path = Path(str(prefix) + ".iqtree")
    if not iqtree_path.exists():
        raise FileNotFoundError(f"Cannot find IQ-TREE report for fitted model: {iqtree_path}")

    text = iqtree_path.read_text(encoding="utf-8", errors="replace")
    for raw in text.splitlines():
        if " -m " not in raw and " -m\t" not in raw:
            continue
        marker = '-m "'
        if marker not in raw:
            continue
        start = raw.index(marker) + len(marker)
        end = raw.find('"', start)
        if end > start:
            fitted = raw[start:end]
            if fitted.startswith("F81+F{") or fitted.startswith("K2P{") or fitted in {"JC", "LG"} or fitted.startswith("GTR{"):
                return fitted

    raise ValueError(f"Could not extract fitted model for {requested_model} from {iqtree_path}")


def write_missing_rows(writer, base, scenario, formula, fitted_evaluation_model):
    row = dict(base)
    row.update(
        {
            "schema_version": SCHEMA_VERSION,
            "scenario": scenario,
            "formula": formula,
            "implementation": "iqtree_native",
            "fitted_evaluation_model": fitted_evaluation_model,
            "target_found": 0,
            "branch_id": "",
            "left_taxa": "",
            "right_taxa": "",
            "valid_sites": "",
            "skipped_sites": "",
            "branch_length_used": "",
            "satC": "",
            "satVar": "",
            "satSE": "",
            "satZ": "",
            "satP": "",
            "decision_unadjusted": "",
            "alpha_taxon_bonf": "",
            "p_taxon_bonf": "",
            "decision_taxon_bonf": "",
            "fdr_by": "",
            "decision_fdr": "",
            "modes": "",
            "eigenvalues": "",
            "weights": "",
        }
    )
    writer.writerow(row)


def write_native_formula_rows(
    writer,
    branch_writer,
    fdr_writer,
    base,
    scenario,
    stat_path,
    fitted_evaluation_model,
    target_taxa,
    alpha,
):
    native_rows = read_native_rows(stat_path)
    target_rows = index_target_rows(native_rows, target_taxa)

    for formula in FORMULAS:
        native = target_rows.get(formula)
        if native is None:
            write_missing_rows(writer, base, scenario, formula, fitted_evaluation_model)
            continue
        row = dict(base)
        row.update(
            {
                "schema_version": SCHEMA_VERSION,
                "scenario": scenario,
                "formula": formula,
                "implementation": "iqtree_native",
                "fitted_evaluation_model": fitted_evaluation_model,
                "target_found": 1,
                "branch_id": native.get("ID", ""),
                "left_taxa": native.get("LeftTaxa", ""),
                "right_taxa": native.get("RightTaxa", ""),
                "valid_sites": native.get("ValidSites", ""),
                "skipped_sites": native.get("SkippedSites", ""),
                "branch_length_used": native.get("Length", ""),
                "alpha": native.get("Alpha", alpha),
                "satC": native.get("satC", ""),
                "satVar": native.get("satVar", ""),
                "satSE": native.get("satSE", ""),
                "satZ": native.get("satZ", ""),
                "satP": native.get("satP", ""),
                "decision_unadjusted": native.get("Decision", ""),
                "alpha_taxon_bonf": native.get("AlphaTaxonBonf", ""),
                "p_taxon_bonf": target_detail(native, alpha)["p_taxon_bonf"],
                "decision_taxon_bonf": native.get("DecisionTaxonBonf", ""),
                "fdr_by": native.get("FDR_BY", ""),
                "decision_fdr": native.get("DecisionFDR", ""),
                "modes": native.get("Modes", ""),
                "eigenvalues": native.get("Eigenvalues", ""),
                "weights": native.get("Weights", ""),
            }
        )
        writer.writerow(row)

    # Persist the complete pooled branch family once per native analysis. This
    # is the auditable denominator for empirical FDR/FDP summaries.
    for formula in FORMULAS:
        formula_rows = [row for row in native_rows if row["Formula"] == formula]
        for branch in branch_audit_detail(formula_rows, target_taxa):
            branch_writer.writerow(
                {
                    "schema_version": SCHEMA_VERSION,
                    "tree_case": base["tree_case"],
                    "simulation_model": base["simulation_model"],
                    "evaluation_model": base["evaluation_model"],
                    "nsites": base["nsites"],
                    "branch_length": base["branch_length"],
                    "replicate": base["replicate"],
                    "seed": base["seed"],
                    "scenario": scenario,
                    "formula": formula,
                    **branch,
                }
            )

    for formula in FORMULAS:
        audit = fdr_family_audit(native_rows, formula)
        fdr_writer.writerow(
            {
                "schema_version": SCHEMA_VERSION,
                "tree_case": base["tree_case"],
                "simulation_model": base["simulation_model"],
                "evaluation_model": base["evaluation_model"],
                "nsites": base["nsites"],
                "branch_length": base["branch_length"],
                "replicate": base["replicate"],
                "seed": base["seed"],
                "scenario": scenario,
                **audit,
            }
        )


def run_case(
    iqtree,
    seqgen,
    indelible,
    simulator,
    pools,
    outdir,
    writer,
    branch_writer,
    fdr_writer,
    timing_writer,
    tree_case,
    simulation_model,
    evaluation_model,
    nsites,
    branch_length,
    rep,
    seed,
    alpha,
    scenario_set,
):
    started = time.perf_counter()
    case_dir = (
        outdir
        / "runs"
        / tree_case
        / f"sim_{sanitize_model(simulation_model)}"
        / f"eval_{sanitize_model(evaluation_model)}"
        / f"n{nsites}"
        / f"b{branch_length:.2f}"
        / f"r{rep:04d}"
    )
    case_dir.mkdir(parents=True, exist_ok=True)
    tree_file = case_dir / "true.tree"
    rng = random.Random(seed)
    write_tree(tree_case, branch_length, tree_file, rng, pools)

    sim_prefix = case_dir / "sim"
    simulation_started = time.perf_counter()
    alignment = simulate_alignment(iqtree, seqgen, indelible, simulator, tree_file, simulation_model, nsites, seed, sim_prefix)
    simulation_seconds = time.perf_counter() - simulation_started
    target_taxa = target_taxa_for_case(tree_case)
    scenarios = set(selected_scenarios(simulation_model, evaluation_model, scenario_set))
    base = {
        "tree_case": tree_case,
        "simulation_model": simulation_model,
        "evaluation_model": evaluation_model,
        "nsites": nsites,
        "branch_length": branch_length,
        "replicate": rep,
        "seed": seed,
        "simulator": simulator,
        "branch_length_source": "evonaps_table" if pools.get("internal") and pools.get("external") else "fixed_0.1_smoke",
        "target_split": ",".join(target_taxa),
        "alpha": alpha,
        "schema_version": SCHEMA_VERSION,
    }
    scenario_seconds = {TRUE_FIXED: 0.0, TRUE_ML_LENGTHS: 0.0, ML_TREE: 0.0}

    if TRUE_FIXED in scenarios:
        analysis_started = time.perf_counter()
        prefix = case_dir / "true_fixed"
        _tree_path, stat_path = run_satute(
            iqtree, alignment, prefix, evaluation_model, tree_file, True, seed
        )
        fitted_model = fitted_model_for_compute(prefix, evaluation_model)
        write_native_formula_rows(
            writer,
            branch_writer,
            fdr_writer,
            base,
            TRUE_FIXED,
            stat_path,
            fitted_model,
            target_taxa,
            alpha,
        )
        scenario_seconds[TRUE_FIXED] = time.perf_counter() - analysis_started

    if TRUE_ML_LENGTHS in scenarios:
        analysis_started = time.perf_counter()
        prefix = case_dir / "true_ml_lengths"
        _tree_path, stat_path = run_satute(
            iqtree, alignment, prefix, evaluation_model, tree_file, False, seed
        )
        fitted_model = fitted_model_for_compute(prefix, evaluation_model)
        write_native_formula_rows(
            writer,
            branch_writer,
            fdr_writer,
            base,
            TRUE_ML_LENGTHS,
            stat_path,
            fitted_model,
            target_taxa,
            alpha,
        )
        scenario_seconds[TRUE_ML_LENGTHS] = time.perf_counter() - analysis_started

    if ML_TREE in scenarios:
        analysis_started = time.perf_counter()
        prefix = case_dir / "ml_tree"
        _tree_path, stat_path = run_satute(
            iqtree, alignment, prefix, evaluation_model, None, False, seed
        )
        fitted_model = fitted_model_for_compute(prefix, evaluation_model)
        write_native_formula_rows(
            writer,
            branch_writer,
            fdr_writer,
            base,
            ML_TREE,
            stat_path,
            fitted_model,
            target_taxa,
            alpha,
        )
        scenario_seconds[ML_TREE] = time.perf_counter() - analysis_started

    timing_writer.writerow(
        {
            "schema_version": SCHEMA_VERSION,
            "tree_case": tree_case,
            "simulation_model": simulation_model,
            "evaluation_model": evaluation_model,
            "nsites": nsites,
            "branch_length": branch_length,
            "replicate": rep,
            "seed": seed,
            "simulation_cache_hit": 0,
            "simulation_seconds": f"{simulation_seconds:.6f}",
            "true_fixed_seconds": f"{scenario_seconds[TRUE_FIXED]:.6f}",
            "true_ml_lengths_seconds": f"{scenario_seconds[TRUE_ML_LENGTHS]:.6f}",
            "ml_tree_seconds": f"{scenario_seconds[ML_TREE]:.6f}",
            "total_seconds": f"{time.perf_counter() - started:.6f}",
        }
    )


def aggregate(detail_path, summary_path):
    # Compatibility wrapper for callers outside the driver.  The implementation
    # lives in the schema-aware module so all consumers use the same rules.
    aggregate_detail(detail_path, summary_path)


def sanitize_model(model):
    return model.replace("{", "_").replace("}", "_").replace("+", "_").replace(",", "_").replace(".", "p")


def task_key_from_row(row):
    return (
        row["tree_case"],
        row["simulation_model"],
        row["evaluation_model"],
        str(row["nsites"]),
        f"{float(row['branch_length']):.10g}",
        str(row["replicate"]),
    )


def completed_tasks(
    detail_path,
    scenario_set,
    fieldnames=None,
    clean_incomplete=False,
    auxiliary_paths=None,
):
    if not detail_path.exists():
        return set()
    rows = []
    seen = defaultdict(set)
    expected = {}
    expected_pairs = {}
    with open(detail_path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fieldnames or BASE_DETAIL_FIELDS):
            raise ValueError(f"Unexpected benchmark detail schema in {detail_path}")
        for row in reader:
            rows.append(row)
            key = task_key_from_row(row)
            pair = (row["scenario"], row["formula"])
            seen[key].add(pair)
            expected[key] = expected_rows_for_task(row["simulation_model"], row["evaluation_model"], scenario_set)
            expected_pairs[key] = {
                (scenario, formula)
                for scenario in selected_scenarios(
                    row["simulation_model"], row["evaluation_model"], scenario_set
                )
                for formula in FORMULAS
            }
    auxiliary_paths = auxiliary_paths or {}
    fdr_seen = defaultdict(set)
    fdr_branch_expected = defaultdict(int)
    branch_counts = defaultdict(int)
    timing_seen = set()
    for path in auxiliary_paths.get("fdr", []):
        if not path.exists():
            continue
        with path.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                key = task_key_from_row(row)
                fdr_seen[key].add((row["scenario"], row["formula"]))
                fdr_branch_expected[key] += int(row["family_size"])
    for path in auxiliary_paths.get("branch", []):
        if not path.exists():
            continue
        with path.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                branch_counts[task_key_from_row(row)] += 1
    for path in auxiliary_paths.get("timing", []):
        if not path.exists():
            continue
        with path.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                timing_seen.add(task_key_from_row(row))

    complete = {
        key for key, pairs in seen.items()
        if len(pairs) == expected.get(key, 0)
        and pairs == expected_pairs.get(key, set())
        and (not auxiliary_paths or fdr_seen[key] == pairs)
        and (not auxiliary_paths or branch_counts[key] == fdr_branch_expected[key])
        and (not auxiliary_paths or key in timing_seen)
    }
    if clean_incomplete and fieldnames is not None:
        def clean_table(path, expected_fields):
            if not path.exists():
                return
            temporary = path.with_name(f".{path.name}.tmp.clean")
            with path.open("r", encoding="utf-8", newline="") as source, temporary.open(
                "w", encoding="utf-8", newline=""
            ) as destination:
                reader = csv.DictReader(source, delimiter="\t")
                writer = csv.DictWriter(destination, fieldnames=expected_fields, delimiter="\t")
                writer.writeheader()
                for row in reader:
                    if task_key_from_row(row) in complete:
                        writer.writerow(row)
            os.replace(temporary, path)

        clean_table(detail_path, fieldnames)
        clean_table(auxiliary_paths.get("branch_path"), BRANCH_AUDIT_FIELDS)
        clean_table(auxiliary_paths.get("fdr_path"), FDR_AUDIT_FIELDS)
        clean_table(auxiliary_paths.get("timing_path"), TIMING_FIELDS)
    return complete


def main():
    parser = argparse.ArgumentParser(
        description="Paired native IQ-TREE dominant versus eigenvalue-weighted SatuTe simulations."
    )
    parser.add_argument("--iqtree", default=str(IQTREE_DEFAULT))
    parser.add_argument("--outdir", default="/tmp/iqtree-satute-head-to-head")
    parser.add_argument("--reps", type=int, default=1)
    parser.add_argument("--site-lengths", default="100")
    parser.add_argument("--branch-lengths", default="0.1,1.0")
    parser.add_argument("--tree-cases", default="five_external")
    parser.add_argument(
        "--simulation-models",
        default="JC",
        help=(
            "Comma-separated simulation model aliases, e.g. JC, GTR_PF06346, "
            "LG, WAG, JTT, Q.PFAM and their *_G4 variants."
        ),
    )
    parser.add_argument(
        "--evaluation-models",
        default="",
        help="Comma-separated evaluation model aliases. Defaults to the simulation models. Use JC,K2P,F81 for misspecification.",
    )
    parser.add_argument(
        "--model-pairs",
        default="",
        help="Exact comma-separated SIM:EVAL model pairs. Use this for matched extension designs such as GTR_SKEW_FREQ:GTR_SKEW_FREQ.",
    )
    parser.add_argument(
        "--scenario-set",
        choices=["fig2", "misspecification", "all"],
        default="fig2",
        help="Scenario contract: fig2 reproduces the main JC simulations; misspecification follows the supplementary GTR setup.",
    )
    parser.add_argument("--paper-grid", action="store_true", help="Use paper branch lengths and site lengths 100,1000,10000.")
    parser.add_argument("--simulator", choices=["alisim", "seq-gen", "indelible"], default="alisim")
    parser.add_argument("--seqgen", default=shutil.which("seq-gen") or shutil.which("seqgen") or "")
    parser.add_argument("--indelible", default=shutil.which("indelible") or shutil.which("INDELible") or "")
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
    if outdir.exists() and not args.resume and any(outdir.iterdir()):
        raise SystemExit(
            f"Refusing to overwrite non-empty output directory {outdir}; "
            "use --resume for a compatible run or choose a new --outdir"
        )
    outdir.mkdir(parents=True, exist_ok=True)
    if args.shard_count < 1:
        raise SystemExit("--shard-count must be >= 1")
    if args.shard_index < 0 or args.shard_index >= args.shard_count:
        raise SystemExit("--shard-index must satisfy 0 <= shard-index < shard-count")

    branch_lengths = parse_csv_numbers(PAPER_BRANCH_LENGTHS if args.paper_grid else args.branch_lengths, float)
    site_lengths = parse_csv_numbers("100,1000,10000" if args.paper_grid else args.site_lengths, int)
    tree_cases = [value.strip() for value in args.tree_cases.split(",") if value.strip()]
    if args.model_pairs:
        try:
            model_pairs = parse_model_pairs(args.model_pairs)
        except ValueError as exc:
            raise SystemExit(str(exc)) from exc
    else:
        simulation_aliases = parse_csv_text(args.simulation_models)
        evaluation_aliases = parse_csv_text(args.evaluation_models) or simulation_aliases
        simulation_models = [resolve_model_alias(alias) for alias in simulation_aliases]
        evaluation_models = [resolve_model_alias(alias) for alias in evaluation_aliases]
        model_pairs = [(simulation_model, evaluation_model) for simulation_model in simulation_models for evaluation_model in evaluation_models]
    pools = load_branch_length_pools(args.evonaps_branch_lengths)

    detail_path = outdir / "head_to_head_detail.tsv"
    summary_path = outdir / "head_to_head_summary.tsv"
    branch_path = outdir / "head_to_head_branch_audit.tsv"
    fdr_path = outdir / "head_to_head_fdr_audit.tsv"
    timing_path = outdir / "head_to_head_timing.tsv"
    fieldnames = list(BASE_DETAIL_FIELDS)
    for path, schema in (
        (detail_path, BASE_DETAIL_FIELDS),
        (summary_path, SUMMARY_FIELDS),
        (branch_path, BRANCH_AUDIT_FIELDS),
        (fdr_path, FDR_AUDIT_FIELDS),
        (timing_path, TIMING_FIELDS),
    ):
        if path is not None and args.resume:
            validate_tsv_header(path, schema)
    done = (
        completed_tasks(
            detail_path,
            args.scenario_set,
            fieldnames,
            clean_incomplete=True,
            auxiliary_paths={
                "branch": [branch_path],
                "fdr": [fdr_path],
                "timing": [timing_path],
                "branch_path": branch_path,
                "fdr_path": fdr_path,
                "timing_path": timing_path,
            },
        )
        if args.resume
        else set()
    )
    modes = {
        "detail": "a" if args.resume and detail_path.exists() and detail_path.stat().st_size > 0 else "w",
        "branch": "a" if args.resume and branch_path.exists() and branch_path.stat().st_size > 0 else "w",
        "fdr": "a" if args.resume and fdr_path.exists() and fdr_path.stat().st_size > 0 else "w",
        "timing": "a" if args.resume and timing_path.exists() and timing_path.stat().st_size > 0 else "w",
    }
    with ExitStack() as stack:
        detail_handle = stack.enter_context(open(detail_path, modes["detail"], encoding="utf-8", newline=""))
        branch_handle = stack.enter_context(open(branch_path, modes["branch"], encoding="utf-8", newline=""))
        fdr_handle = stack.enter_context(open(fdr_path, modes["fdr"], encoding="utf-8", newline=""))
        timing_handle = stack.enter_context(open(timing_path, modes["timing"], encoding="utf-8", newline=""))
        writer = csv.DictWriter(detail_handle, fieldnames=fieldnames, delimiter="\t")
        branch_writer = csv.DictWriter(branch_handle, fieldnames=BRANCH_AUDIT_FIELDS, delimiter="\t")
        fdr_writer = csv.DictWriter(fdr_handle, fieldnames=FDR_AUDIT_FIELDS, delimiter="\t")
        timing_writer = csv.DictWriter(timing_handle, fieldnames=TIMING_FIELDS, delimiter="\t")
        if modes["detail"] == "w":
            writer.writeheader()
        if modes["branch"] == "w":
            branch_writer.writeheader()
        if modes["fdr"] == "w":
            fdr_writer.writeheader()
        if modes["timing"] == "w":
            timing_writer.writeheader()
        task_index = 0
        selected_tasks = 0
        skipped_tasks = 0
        for tree_case in tree_cases:
            for simulation_model, evaluation_model in model_pairs:
                for nsites in site_lengths:
                    for branch_length in branch_lengths:
                        for rep in range(1, args.reps + 1):
                            current = task_index
                            task_index += 1
                            if current % args.shard_count != args.shard_index:
                                continue
                            selected_tasks += 1
                            task_key = (
                                tree_case,
                                simulation_model,
                                evaluation_model,
                                str(nsites),
                                f"{branch_length:.10g}",
                                str(rep),
                            )
                            if task_key in done:
                                skipped_tasks += 1
                                continue
                            seed = simulation_seed(900000, tree_case, simulation_model, nsites, branch_length, rep)
                            run_case(
                                iqtree,
                                args.seqgen,
                                args.indelible,
                                args.simulator,
                                pools,
                                outdir,
                                writer,
                                branch_writer,
                                fdr_writer,
                                timing_writer,
                                tree_case,
                                simulation_model,
                                evaluation_model,
                                nsites,
                                branch_length,
                                rep,
                                seed,
                                args.alpha,
                                args.scenario_set,
                            )
                            detail_handle.flush()
                            branch_handle.flush()
                            fdr_handle.flush()
                            timing_handle.flush()

    aggregate_detail(detail_path, summary_path)
    print(f"Shard:   {args.shard_index}/{args.shard_count}")
    print(f"Tasks:   {selected_tasks}/{task_index}")
    if args.resume:
        print(f"Skipped: {skipped_tasks}")
    print(f"Detail:  {detail_path}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
