#!/usr/bin/env python3

import argparse
import math
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


SATUTE_SRC = Path(__file__).resolve().parents[2] / "src"
if str(SATUTE_SRC) not in sys.path:
    sys.path.insert(0, str(SATUTE_SRC))

from satute_analysis.multiple_testing import annotate_satute_rows


STATE_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}
FORMULA_VARIANTS = ("dominant", "eigenvalue_weighted")


@dataclass
class TreeNode:
    name: str | None = None
    length: float = 0.0
    children: list["TreeNode"] = field(default_factory=list)
    parent: "TreeNode | None" = None
    node_id: int = -1

    def is_leaf(self):
        return not self.children


@dataclass
class SplitEdge:
    left_node: TreeNode
    left_dad: TreeNode
    right_node: TreeNode
    right_dad: TreeNode
    branch_length: float
    left_taxa: tuple[str, ...]
    right_taxa: tuple[str, ...]


def parse_fasta(path):
    sequences = {}
    name = None
    chunks = []
    with Path(path).open("r", encoding="utf-8") as handle:
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
    if len(sequences) < 2:
        raise ValueError("Reference calculation requires at least two taxa")
    lengths = {len(seq) for seq in sequences.values()}
    if len(lengths) != 1:
        raise ValueError("All sequences must have the same length")
    return sequences


def split_label(taxa):
    return ",".join(sorted(taxa))


def parse_split(value, taxa):
    left = tuple(item.strip() for item in value.split(",") if item.strip())
    if not left:
        raise ValueError("Split must name at least one taxon")
    if len(set(left)) != len(left):
        raise ValueError(f"Duplicate taxon in split: {value}")
    unknown = sorted(set(left) - set(taxa))
    if unknown:
        raise ValueError(f"Split contains unknown taxa: {unknown}")
    if len(left) >= len(taxa):
        raise ValueError("Split must leave at least one taxon on the opposite side")
    right = tuple(sorted(set(taxa) - set(left)))
    return tuple(sorted(left)), right


def strip_newick_comments(text):
    out = []
    depth = 0
    for char in text:
        if char == "[":
            depth += 1
            continue
        if char == "]" and depth > 0:
            depth -= 1
            continue
        if depth == 0:
            out.append(char)
    return "".join(out)


def parse_newick(path):
    text = strip_newick_comments(Path(path).read_text(encoding="utf-8")).strip()
    if not text.endswith(";"):
        raise ValueError(f"Newick tree must end with ';': {path}")
    text = "".join(text.split())
    index = 0

    def parse_name():
        nonlocal index
        start = index
        while index < len(text) and text[index] not in "(),:;":
            index += 1
        return text[start:index] or None

    def parse_length():
        nonlocal index
        if index >= len(text) or text[index] != ":":
            return 0.0
        index += 1
        start = index
        while index < len(text) and text[index] not in ",);":
            index += 1
        try:
            return float(text[start:index])
        except ValueError as exc:
            raise ValueError(f"Invalid branch length in {path}: {text[start:index]}") from exc

    def parse_subtree():
        nonlocal index
        if index >= len(text):
            raise ValueError(f"Unexpected end of Newick tree in {path}")

        if text[index] == "(":
            index += 1
            node = TreeNode()
            while True:
                child = parse_subtree()
                child.parent = node
                node.children.append(child)
                if index >= len(text):
                    raise ValueError(f"Unclosed internal node in {path}")
                if text[index] == ",":
                    index += 1
                    continue
                if text[index] == ")":
                    index += 1
                    break
                raise ValueError(f"Unexpected character in {path}: {text[index]}")
            node.name = parse_name()
            node.length = parse_length()
            return node

        name = parse_name()
        if not name:
            raise ValueError(f"Expected leaf name in {path}")
        node = TreeNode(name=name)
        node.length = parse_length()
        return node

    root = parse_subtree()
    if index >= len(text) or text[index] != ";":
        raise ValueError(f"Unexpected trailing Newick content in {path}: {text[index:]}")
    assign_node_ids(root)
    return root


def assign_node_ids(root):
    counter = 0
    stack = [root]
    while stack:
        node = stack.pop()
        node.node_id = counter
        counter += 1
        stack.extend(reversed(node.children))


def iter_neighbors(node):
    if node.parent is not None:
        yield node.parent, node.length
    for child in node.children:
        yield child, child.length


def tree_edges(root):
    edges = []
    stack = [root]
    while stack:
        node = stack.pop()
        for child in node.children:
            edges.append((node, child, child.length))
            stack.append(child)
    return edges


def collect_tree_taxa(root):
    taxa = []
    stack = [root]
    while stack:
        node = stack.pop()
        if node.is_leaf():
            if not node.name:
                raise ValueError("Encountered unnamed leaf in tree")
            taxa.append(node.name)
        else:
            stack.extend(reversed(node.children))
    if len(taxa) != len(set(taxa)):
        duplicated = sorted(name for name in set(taxa) if taxa.count(name) > 1)
        raise ValueError(f"Duplicated taxon names in tree: {duplicated}")
    return tuple(sorted(taxa))


def component_taxa(node, blocked):
    taxa = []

    def visit(current, parent):
        if current.is_leaf():
            taxa.append(current.name)
            return
        for neighbor, _length in iter_neighbors(current):
            if neighbor is parent:
                continue
            visit(neighbor, current)

    visit(node, blocked)
    return frozenset(taxa)


def find_split_edge(root, requested_taxa):
    all_taxa = frozenset(collect_tree_taxa(root))
    requested = frozenset(requested_taxa)
    for parent, child, length in tree_edges(root):
        child_side = component_taxa(child, parent)
        parent_side = all_taxa - child_side
        if child_side == requested:
            return SplitEdge(
                left_node=child,
                left_dad=parent,
                right_node=parent,
                right_dad=child,
                branch_length=length,
                left_taxa=tuple(sorted(child_side)),
                right_taxa=tuple(sorted(parent_side)),
            )
        if parent_side == requested:
            return SplitEdge(
                left_node=parent,
                left_dad=child,
                right_node=child,
                right_dad=parent,
                branch_length=length,
                left_taxa=tuple(sorted(parent_side)),
                right_taxa=tuple(sorted(child_side)),
            )
    raise ValueError(f"Requested split {split_label(requested)} is not a branch in the tree")


def all_split_edges(root):
    """Return every unique edge of an explicitly unrooted Newick tree.

    A two-child top-level node represents a rooted or degree-two-rooted tree and
    splits one unrooted edge into two lengths.  Rejecting that representation
    avoids silently using only half of the focal branch length.
    """

    if len(root.children) == 2:
        raise ValueError(
            "--all-branches requires an unrooted Newick representation with a "
            "top-level trifurcation, not a two-child root"
        )

    all_taxa = frozenset(collect_tree_taxa(root))
    split_edges = []
    seen = set()
    for parent, child, length in tree_edges(root):
        child_side = component_taxa(child, parent)
        parent_side = all_taxa - child_side
        if len(child_side) <= len(parent_side):
            left_node, left_dad = child, parent
            right_node, right_dad = parent, child
            left_side, right_side = child_side, parent_side
        else:
            left_node, left_dad = parent, child
            right_node, right_dad = child, parent
            left_side, right_side = parent_side, child_side

        key = (tuple(sorted(left_side)), tuple(sorted(right_side)))
        reverse_key = (key[1], key[0])
        canonical_key = min(key, reverse_key)
        if canonical_key in seen:
            continue
        seen.add(canonical_key)
        split_edges.append(
            SplitEdge(
                left_node=left_node,
                left_dad=left_dad,
                right_node=right_node,
                right_dad=right_dad,
                branch_length=length,
                left_taxa=tuple(sorted(left_side)),
                right_taxa=tuple(sorted(right_side)),
            )
        )

    return sorted(split_edges, key=lambda edge: (len(edge.left_taxa), edge.left_taxa))


def parse_model(model):
    if model.startswith("JC"):
        return np.ones(6, dtype=float), np.full(4, 0.25, dtype=float)

    match = re.search(r"GTR\{([^}]+)\}\+F\{([^}]+)\}", model)
    if not match:
        raise ValueError(f"Unsupported model syntax: {model}")
    rates = np.array([float(value) for value in match.group(1).split(",")], dtype=float)
    frequencies = np.array([float(value) for value in match.group(2).split(",")], dtype=float)
    frequencies /= np.sum(frequencies)
    if len(rates) != 6 or len(frequencies) != 4:
        raise ValueError(f"Unsupported model dimensions: {model}")
    return rates, frequencies


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
    transition = eigenvectors @ np.diag(np.exp(evals * length)) @ inv_eigenvectors
    close_negative = (transition < 0.0) & (transition > -1e-14)
    transition[close_negative] = 0.0
    return transition


def formula_modes_and_weights(evals, formula, effective_branch_length):
    zero_index = int(np.argmin(np.abs(evals)))
    nonzero = [i for i, value in enumerate(evals) if i != zero_index and abs(value) > 1e-10]
    if not nonzero:
        raise ValueError("No non-stationary eigenmodes found")

    dominant = max(evals[i] for i in nonzero)
    if formula == "dominant":
        tol = max(1e-8, abs(dominant) * 1e-6)
        modes = [i for i in nonzero if abs(evals[i] - dominant) <= tol]
        weights = np.ones(len(modes), dtype=float)
    elif formula == "all_unweighted":
        modes = nonzero
        weights = np.ones(len(modes), dtype=float)
    elif formula == "eigenvalue_weighted":
        modes = nonzero
        weights = np.exp((evals[modes] - dominant) * effective_branch_length)
    else:
        raise ValueError(f"Unsupported formula: {formula}")
    return modes, weights


def side_partial(sequences, site, node, dad, transition_builder):
    if node.is_leaf():
        state = STATE_INDEX.get(sequences[node.name][site])
        if state is None:
            return np.ones(4, dtype=float)
        partial = np.zeros(4, dtype=float)
        partial[state] = 1.0
        return partial

    partial = np.ones(4, dtype=float)
    for child, length in iter_neighbors(node):
        if child is dad:
            continue
        child_partial = side_partial(sequences, site, child, node, transition_builder)
        transition = transition_builder(length)
        partial *= transition @ child_partial
    return partial


def posterior(partial, pi):
    values = partial * pi
    total = float(np.sum(values))
    if not math.isfinite(total) or total <= 0.0:
        return None
    return values / total


def compute_reference(sequences, split_edge, sites, model, formula, rate_multiplier, alpha):
    q, pi = build_q(model)
    evals, eigenvectors, inv_eigenvectors = reversible_eigendecomposition(q, pi)
    effective_branch_length = split_edge.branch_length * rate_multiplier
    modes, weights = formula_modes_and_weights(evals, formula, effective_branch_length)
    transition_cache = {}

    def build_transition(length):
        scaled_length = length * rate_multiplier
        if scaled_length not in transition_cache:
            transition_cache[scaled_length] = transition_matrix(
                evals,
                eigenvectors,
                inv_eigenvectors,
                scaled_length,
            )
        return transition_cache[scaled_length]

    coherence_sum = 0.0
    left_second = np.zeros((len(modes), len(modes)), dtype=float)
    right_second = np.zeros((len(modes), len(modes)), dtype=float)
    valid_sites = 0

    for site in sites:
        left_post = posterior(
            side_partial(sequences, site, split_edge.left_node, split_edge.left_dad, build_transition),
            pi,
        )
        right_post = posterior(
            side_partial(sequences, site, split_edge.right_node, split_edge.right_dad, build_transition),
            pi,
        )
        if left_post is None or right_post is None:
            continue

        left_factors = np.array([float(np.dot(eigenvectors[:, mode], left_post)) for mode in modes])
        right_factors = np.array([float(np.dot(eigenvectors[:, mode], right_post)) for mode in modes])
        coherence_sum += float(np.dot(weights, left_factors * right_factors))
        left_second += np.outer(left_factors, left_factors)
        right_second += np.outer(right_factors, right_factors)
        valid_sites += 1

    if valid_sites == 0:
        raise ValueError("No valid sites for the requested reference calculation")

    sat_c = coherence_sum / valid_sites
    sat_var = float(
        np.sum(np.outer(weights, weights) * (left_second / valid_sites) * (right_second / valid_sites))
    )
    sat_se = math.sqrt(sat_var / valid_sites)
    sat_z = sat_c / sat_se
    sat_p = 0.5 * math.erfc(sat_z / math.sqrt(2.0))
    information_fraction = spectral_information_fraction(
        evals,
        split_edge.branch_length,
        [{"rate": rate_multiplier, "proportion": 1.0}],
    )
    return {
        "valid_sites": valid_sites,
        "satC": sat_c,
        "satVar": sat_var,
        "satSE": sat_se,
        "satZ": sat_z,
        "satP": sat_p,
        "InformationFraction": information_fraction,
        "SaturationIndex": 1.0 - information_fraction,
        "Decision": "informative" if sat_p <= alpha else "saturated",
        "Modes": ",".join(str(mode) for mode in modes),
        "Eigenvalues": ",".join(format_float(float(evals[mode])) for mode in modes),
        "Weights": ",".join(format_float(float(weight)) for weight in weights),
    }


def parse_rate_file(path, nsites):
    if path is None:
        return {"pooled": {"rate": 1.0, "sites": list(range(nsites))}}

    header = None
    categories = defaultdict(lambda: {"rate": None, "sites": []})
    with Path(path).open("r", encoding="utf-8") as handle:
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
            site = int(row["Site"]) - 1
            if site < 0 or site >= nsites:
                raise ValueError(f"Site index outside alignment length in {path}: {site + 1}")
            cat = row["Cat"]
            categories[cat]["rate"] = float(row["C_Rate"])
            categories[cat]["sites"].append(site)
    return dict(categories)


def spectral_information_fraction(evals, branch_length, rate_categories):
    zero_index = int(np.argmin(np.abs(evals)))
    modes = [i for i, value in enumerate(evals) if i != zero_index and abs(value) > 1e-10]
    if not modes:
        raise ValueError("No nonstationary eigenmodes for the saturation scale")
    total_proportion = sum(category["proportion"] for category in rate_categories)
    if total_proportion <= 0.0:
        raise ValueError("Rate-category proportions do not sum to a positive value")
    information = 0.0
    for category in rate_categories:
        effective_length = branch_length * category["rate"]
        category_information = float(np.mean(np.exp(2.0 * evals[modes] * effective_length)))
        information += category["proportion"] * category_information
    return information / total_proportion


def parse_iqtree_rate_categories(path):
    categories = []
    in_table = False
    with Path(path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line == "Category  Relative_rate  Proportion":
                in_table = True
                continue
            if not in_table:
                continue
            fields = line.split()
            if len(fields) != 3 or not fields[0].isdigit():
                if categories:
                    break
                continue
            categories.append(
                {
                    "label": fields[0],
                    "rate": float(fields[1]),
                    "proportion": float(fields[2]),
                }
            )
    if not categories:
        raise ValueError(f"No rate-category table found in {path}")
    total = sum(category["proportion"] for category in categories)
    if total <= 0.0:
        raise ValueError(f"Rate-category proportions do not sum to a positive value in {path}")
    for category in categories:
        category["proportion"] /= total
    return categories


def pool_results(category_results, alpha):
    valid_total = sum(result["valid_sites"] for result in category_results)
    assigned_total = sum(result["assigned_sites"] for result in category_results)
    sat_c = 0.0
    sat_var = 0.0
    for result in category_results:
        weight = result["valid_sites"] / valid_total
        sat_c += weight * result["satC"]
        sat_var += weight * result["satVar"]
    sat_se = math.sqrt(sat_var / valid_total)
    sat_z = sat_c / sat_se
    sat_p = 0.5 * math.erfc(sat_z / math.sqrt(2.0))
    return {
        "assigned_sites": assigned_total,
        "valid_sites": valid_total,
        "satC": sat_c,
        "satVar": sat_var,
        "satSE": sat_se,
        "satZ": sat_z,
        "satP": sat_p,
        "Decision": "informative" if sat_p <= alpha else "saturated",
        "Modes": category_results[0]["Modes"],
        "Eigenvalues": category_results[0]["Eigenvalues"],
        "Weights": "NA",
    }


def compute_split_rows(
    sequences,
    split_edge,
    categories,
    rate_categories_for_scale,
    model,
    evals_for_scale,
    alpha,
):
    pooled_information_fraction = spectral_information_fraction(
        evals_for_scale,
        split_edge.branch_length,
        rate_categories_for_scale,
    )
    split_metadata = {
        "Split": split_label(split_edge.left_taxa),
        "OppositeSplit": split_label(split_edge.right_taxa),
        "LeftTaxa": len(split_edge.left_taxa),
        "RightTaxa": len(split_edge.right_taxa),
    }
    rows = []

    for formula in FORMULA_VARIANTS:
        category_results = []
        category_rows = []
        for category, category_info in sorted(categories.items(), key=lambda item: item[0]):
            reference = compute_reference(
                sequences,
                split_edge,
                category_info["sites"],
                model,
                formula,
                category_info["rate"],
                alpha,
            )
            reference = {
                **reference,
                "assigned_sites": len(category_info["sites"]),
            }
            category_results.append(reference)
            if category != "pooled":
                category_rows.append(
                    {
                        **split_metadata,
                        "Formula": formula,
                        "RateCategory": category,
                        "RateMultiplier": format_float(category_info["rate"]),
                        "Sites": reference["assigned_sites"],
                        **reference,
                    }
                )

        pooled = (
            category_results[0]
            if list(categories) == ["pooled"]
            else pool_results(category_results, alpha)
        )
        pooled["InformationFraction"] = pooled_information_fraction
        pooled["SaturationIndex"] = 1.0 - pooled_information_fraction
        rows.append(
            {
                **split_metadata,
                "Formula": formula,
                "RateCategory": "pooled",
                "RateMultiplier": "NA",
                "Sites": pooled.get("assigned_sites", pooled["valid_sites"]),
                **pooled,
            }
        )
        rows.extend(category_rows)

    return rows


def format_float(value):
    return f"{value:.10g}"


def format_optional_float(value):
    return "NA" if value is None else format_float(value)


def write_rows(path, rows):
    header = [
        "Split",
        "OppositeSplit",
        "Formula",
        "RateCategory",
        "RateMultiplier",
        "Sites",
        "LeftTaxa",
        "RightTaxa",
        "satC",
        "satVar",
        "satSE",
        "satZ",
        "satP",
        "Alpha",
        "AlphaTaxonBonf",
        "PTaxonBonf",
        "DecisionTaxonBonf",
        "FDR_BY",
        "DecisionFDR",
        "InformationFraction",
        "SaturationIndex",
        "Decision",
        "Modes",
        "Eigenvalues",
        "Weights",
    ]
    lines = ["\t".join(header)]
    for row in rows:
        lines.append(
            "\t".join(
                [
                    row["Split"],
                    row["OppositeSplit"],
                    row["Formula"],
                    row["RateCategory"],
                    row["RateMultiplier"],
                    str(row["Sites"]),
                    str(row["LeftTaxa"]),
                    str(row["RightTaxa"]),
                    format_float(row["satC"]),
                    format_float(row["satVar"]),
                    format_float(row["satSE"]),
                    format_float(row["satZ"]),
                    format_float(row["satP"]),
                    format_float(row["Alpha"]),
                    format_float(row["AlphaTaxonBonf"]),
                    format_optional_float(row["PTaxonBonf"]),
                    row["DecisionTaxonBonf"],
                    format_optional_float(row["FDR_BY"]),
                    row["DecisionFDR"],
                    format_float(row["InformationFraction"]),
                    format_float(row["SaturationIndex"]),
                    row["Decision"],
                    row["Modes"],
                    row["Eigenvalues"],
                    row["Weights"],
                ]
            )
        )
    text = "\n".join(lines) + "\n"
    if path is None:
        print(text, end="")
    else:
        Path(path).write_text(text, encoding="utf-8")


def native_split_for_branch_id(path, branch_id):
    header = None
    branch_id = str(branch_id)
    with Path(path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[0] == "ID":
                header = fields
                continue
            if header is None or not fields[0].isdigit():
                continue
            row = dict(zip(header, fields))
            if row.get("ID") == branch_id and row.get("Split"):
                return row["Split"]
    raise ValueError(f"Branch ID {branch_id} was not found in {path}")


def parse_native_sat_stat(path, split_labels=None, branch_id=None):
    header = None
    rows = {}
    branch_id = None if branch_id is None else str(branch_id)
    with Path(path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[0] == "ID":
                header = fields
                continue
            if header is None or not fields[0].isdigit():
                continue
            row = dict(zip(header, fields))
            if branch_id is not None:
                if row.get("ID") != branch_id:
                    continue
            elif split_labels is not None and row.get("Split") not in split_labels:
                continue
            if split_labels is None:
                rows[(row["Split"], row["Formula"], row["RateCategory"])] = row
            elif row.get("Split") in split_labels:
                rows[(row["Formula"], row["RateCategory"])] = row
    return rows


def mode_count(value):
    if not value:
        return 0
    return len([item for item in value.split(",") if item])


def compare_native_rows(native_path, reference_rows, tolerance, split_labels, branch_id=None):
    native_rows = parse_native_sat_stat(native_path, split_labels, branch_id)
    checked = 0

    for reference in reference_rows:
        key = (reference["Formula"], reference["RateCategory"])
        native = native_rows.get(key)
        if native is None:
            raise AssertionError(f"Missing native SatuTe row for {key}")

        native_sites = int(native["RateSites"])
        if native_sites != reference["Sites"]:
            raise AssertionError(
                f"{key}: native sites {native_sites} differ from reference {reference['Sites']}"
            )

        for field in (
            "satC",
            "satVar",
            "satSE",
            "satZ",
            "satP",
            "InformationFraction",
            "SaturationIndex",
        ):
            observed = float(native[field])
            expected = float(reference[field])
            if not math.isclose(observed, expected, rel_tol=0.0, abs_tol=tolerance):
                raise AssertionError(
                    f"{key} {field}: native {observed} differs from reference {expected}"
                )

        if native["Decision"] != reference["Decision"]:
            raise AssertionError(
                f"{key}: native decision {native['Decision']} differs from "
                f"reference {reference['Decision']}"
            )

        if mode_count(native["Modes"]) != mode_count(reference["Modes"]):
            raise AssertionError(
                f"{key}: native mode count {mode_count(native['Modes'])} differs from "
                f"reference {mode_count(reference['Modes'])}"
            )

        checked += 1

    return checked


def compare_native_all_rows(native_path, reference_rows, tolerance):
    native_rows = parse_native_sat_stat(native_path)
    checked = 0

    for reference in reference_rows:
        key = (reference["Split"], reference["Formula"], reference["RateCategory"])
        native = native_rows.get(key)
        if native is None:
            opposite_key = (
                reference["OppositeSplit"],
                reference["Formula"],
                reference["RateCategory"],
            )
            native = native_rows.get(opposite_key)
        if native is None:
            raise AssertionError(f"Missing native SatuTe row for {key}")

        if int(native["RateSites"]) != reference["Sites"]:
            raise AssertionError(
                f"{key}: native sites {native['RateSites']} differ from reference {reference['Sites']}"
            )

        for field in (
            "satC",
            "satVar",
            "satSE",
            "satZ",
            "satP",
            "AlphaTaxonBonf",
            "InformationFraction",
            "SaturationIndex",
        ):
            observed = float(native[field])
            expected = float(reference[field])
            if not math.isclose(observed, expected, rel_tol=0.0, abs_tol=tolerance):
                raise AssertionError(
                    f"{key} {field}: native {observed} differs from reference {expected}"
                )

        for field in ("Decision", "DecisionTaxonBonf"):
            if native[field] != reference[field]:
                raise AssertionError(
                    f"{key} {field}: native {native[field]} differs from reference {reference[field]}"
                )

        reference_fdr = reference["FDR_BY"]
        if reference_fdr is None:
            if native["FDR_BY"] != "NA" or native["DecisionFDR"] != "not_tested":
                raise AssertionError(f"{key}: native category row entered an FDR family")
        else:
            observed_fdr = float(native["FDR_BY"])
            if not math.isclose(observed_fdr, reference_fdr, rel_tol=0.0, abs_tol=tolerance):
                raise AssertionError(
                    f"{key} FDR_BY: native {observed_fdr} differs from reference {reference_fdr}"
                )
            if native["DecisionFDR"] != reference["DecisionFDR"]:
                raise AssertionError(
                    f"{key} DecisionFDR: native {native['DecisionFDR']} differs from "
                    f"reference {reference['DecisionFDR']}"
                )

        checked += 1

    if checked != len(native_rows):
        raise AssertionError(
            f"Python reference covered {checked} rows but native output contains {len(native_rows)} rows"
        )
    return checked


def main():
    parser = argparse.ArgumentParser(
        description="Independent nucleotide Python reference for phase-1 SatuTe formulas."
    )
    parser.add_argument("--alignment", required=True, help="FASTA alignment with nucleotide taxa")
    parser.add_argument("--tree", required=True, help="Newick tree containing the requested split")
    parser.add_argument("--model", required=True, help="JC or GTR{...}+F{...} model string")
    parser.add_argument(
        "--split",
        help="Comma-separated taxa on one side of one target branch (default: A,B)",
    )
    parser.add_argument(
        "--all-branches",
        action="store_true",
        help="Evaluate every unrooted branch and form tree-wide BY families",
    )
    parser.add_argument("--rate-file", help="Optional IQ-TREE .rate file for category rows")
    parser.add_argument("--iqtree-report", help="IQ-TREE .iqtree report containing category rates and proportions")
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--out", help="Optional output TSV path; stdout is used by default")
    parser.add_argument("--compare-sat-stat", help="Optional native .sat.stat file to compare against")
    parser.add_argument("--branch-id", help="Optional native branch ID to target; requires --compare-sat-stat")
    parser.add_argument("--tolerance", type=float, default=2e-4)
    args = parser.parse_args()

    if not (0.0 < args.alpha < 1.0):
        raise SystemExit("--alpha must be between 0 and 1")
    if args.all_branches and args.split is not None:
        raise SystemExit("--all-branches and --split are mutually exclusive")
    if args.all_branches and args.branch_id is not None:
        raise SystemExit("--all-branches and --branch-id are mutually exclusive")

    sequences = parse_fasta(args.alignment)
    root = parse_newick(args.tree)
    tree_taxa = set(collect_tree_taxa(root))
    sequence_taxa = set(sequences)
    if tree_taxa != sequence_taxa:
        raise SystemExit(
            f"Tree/alignment taxa differ: tree-only={sorted(tree_taxa - sequence_taxa)} "
            f"alignment-only={sorted(sequence_taxa - tree_taxa)}"
        )
    if args.branch_id is not None and not args.compare_sat_stat:
        raise SystemExit("--branch-id requires --compare-sat-stat")

    split_labels = None
    if args.all_branches:
        try:
            split_edges = all_split_edges(root)
        except ValueError as exc:
            raise SystemExit(str(exc)) from exc
    else:
        split_value = args.split or "A,B"
        if args.branch_id is not None:
            try:
                split_value = native_split_for_branch_id(args.compare_sat_stat, args.branch_id)
            except ValueError as exc:
                raise SystemExit(str(exc)) from exc
        try:
            left, right = parse_split(split_value, sequences)
            split_edges = [find_split_edge(root, left)]
        except ValueError as exc:
            raise SystemExit(str(exc)) from exc
        split_labels = {split_label(left), split_label(right)}

    nsites = len(next(iter(sequences.values())))
    categories = parse_rate_file(args.rate_file, nsites)
    if args.rate_file is None:
        rate_categories_for_scale = [{"label": "pooled", "rate": 1.0, "proportion": 1.0}]
    else:
        if not args.iqtree_report:
            raise SystemExit("--iqtree-report is required with --rate-file for the pooled information scale")
        rate_categories_for_scale = parse_iqtree_rate_categories(args.iqtree_report)
    q_for_scale, _pi_for_scale = build_q(args.model)
    evals_for_scale, _eigenvectors_for_scale, _inverse_for_scale = reversible_eigendecomposition(
        q_for_scale,
        _pi_for_scale,
    )
    rows = []
    for split_edge in split_edges:
        rows.extend(
            compute_split_rows(
                sequences,
                split_edge,
                categories,
                rate_categories_for_scale,
                args.model,
                evals_for_scale,
                args.alpha,
            )
        )

    rows = annotate_satute_rows(rows, args.alpha, FORMULA_VARIANTS)
    write_rows(args.out, rows)

    if args.compare_sat_stat:
        if args.all_branches:
            checked = compare_native_all_rows(args.compare_sat_stat, rows, args.tolerance)
        else:
            checked = compare_native_rows(
                args.compare_sat_stat,
                rows,
                args.tolerance,
                split_labels,
                args.branch_id,
            )
        print(
            f"Reference matched {checked} native .sat.stat rows within {args.tolerance:g}",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
