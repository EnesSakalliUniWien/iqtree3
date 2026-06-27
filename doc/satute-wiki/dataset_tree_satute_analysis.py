#!/usr/bin/env python3

import argparse
import csv
import json
import math
import re
import statistics
import subprocess
from collections import defaultdict
from pathlib import Path


IQTREE_FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"


class NewickNode:
    next_id = 0

    def __init__(self, label=""):
        self.id = NewickNode.next_id
        NewickNode.next_id += 1
        self.label = label
        self.children = []
        self.length = None
        self.parent = None


def run_command(cmd, log_path):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w", encoding="utf-8") as log:
        log.write("$ " + " ".join(str(part) for part in cmd) + "\n\n")
        log.flush()
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)


def parse_label(text, idx):
    start = idx
    quote = None
    while idx < len(text):
        ch = text[idx]
        if quote:
            if ch == quote:
                quote = None
            idx += 1
            continue
        if ch in "'\"":
            quote = ch
            idx += 1
            continue
        if ch in ":,();":
            break
        idx += 1
    return text[start:idx].strip().strip("'\""), idx


def parse_length(text, idx):
    while idx < len(text) and text[idx].isspace():
        idx += 1
    if idx >= len(text) or text[idx] != ":":
        return None, idx
    idx += 1
    start = idx
    while idx < len(text) and text[idx] not in ",();":
        idx += 1
    raw = text[start:idx].strip()
    return float(raw) if raw else None, idx


def parse_newick(path):
    text = Path(path).read_text(encoding="utf-8").strip()
    if not text:
        raise ValueError(f"Empty Newick file: {path}")
    NewickNode.next_id = 0
    idx = 0

    def skip_ws():
        nonlocal idx
        while idx < len(text) and text[idx].isspace():
            idx += 1

    def parse_subtree():
        nonlocal idx
        skip_ws()
        if idx >= len(text):
            raise ValueError("Unexpected end of Newick")
        node = NewickNode()
        if text[idx] == "(":
            idx += 1
            while True:
                child = parse_subtree()
                child.parent = node
                node.children.append(child)
                skip_ws()
                if idx < len(text) and text[idx] == ",":
                    idx += 1
                    continue
                if idx < len(text) and text[idx] == ")":
                    idx += 1
                    break
                raise ValueError(f"Unexpected Newick text near {text[idx:idx + 40]!r}")
            node.label, idx = parse_label(text, idx)
        else:
            node.label, idx = parse_label(text, idx)
            if not node.label:
                raise ValueError(f"Missing leaf label near {text[idx:idx + 40]!r}")
        node.length, idx = parse_length(text, idx)
        return node

    root = parse_subtree()
    skip_ws()
    if idx < len(text) and text[idx] == ";":
        idx += 1
    skip_ws()
    if idx != len(text):
        raise ValueError(f"Unexpected trailing Newick text: {text[idx:]!r}")
    return root


def leaf_names(node):
    if not node.children:
        return {node.label}
    names = set()
    for child in node.children:
        names.update(leaf_names(child))
    return names


def iter_nodes(root):
    stack = [root]
    while stack:
        node = stack.pop()
        yield node
        stack.extend(reversed(node.children))


def canonical_split(side, all_taxa):
    side = frozenset(t for t in side if t)
    other = frozenset(all_taxa - side)
    if not side or not other:
        return None
    return tuple(sorted(side if len(side) <= len(other) else other))


def parse_support_label(label, support_mode):
    if not label:
        return "", "", ""
    parts = [part.strip() for part in str(label).split("/") if part.strip()]
    numeric = []
    for part in parts:
        try:
            numeric.append(float(part))
        except ValueError:
            numeric.append(math.nan)
    if len(numeric) >= 2:
        return label, numeric[0], numeric[1]
    if len(numeric) == 1:
        if support_mode == "ufboot":
            return label, "", numeric[0]
        return label, numeric[0], ""
    return label, "", ""


def tree_support_by_split(tree_path, support_mode):
    root = parse_newick(tree_path)
    taxa = leaf_names(root)
    support = {}
    branch_lengths = {}
    for node in iter_nodes(root):
        if node is root:
            continue
        split = canonical_split(leaf_names(node), taxa)
        if split is None:
            continue
        label, sh_alrt, ufboot = parse_support_label(node.label if node.children else "", support_mode)
        support[split] = {
            "support_label": label,
            "sh_alrt": sh_alrt,
            "ufboot": ufboot,
        }
        branch_lengths[split] = node.length
    return support, branch_lengths, taxa


def split_from_sat_row(row, all_taxa):
    split = row.get("Split", "")
    if not split:
        return None
    side = {part.strip() for part in split.split(",") if part.strip()}
    return canonical_split(side, all_taxa)


def parse_sat_stat(path):
    rows = []
    header = None
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[0] == "ID":
                header = fields
                continue
            if header is None:
                continue
            row = dict(zip(header, fields))
            if row.get("RateCategory", "pooled") != "pooled":
                continue
            rows.append(row)
    return rows


def extract_metric(text, label):
    match = re.search(rf"{re.escape(label)}:\s*({IQTREE_FLOAT})", text)
    return match.group(1) if match else ""


def parse_iqtree_report(path):
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    return {
        "log_likelihood": extract_metric(text, "Log-likelihood of the tree"),
        "unconstrained_log_likelihood": extract_metric(text, "Unconstrained log-likelihood (without tree)"),
        "aic": extract_metric(text, "Akaike information criterion (AIC) score"),
        "aicc": extract_metric(text, "Corrected Akaike information criterion (AICc) score"),
        "bic": extract_metric(text, "Bayesian information criterion (BIC) score"),
    }


def build_iqtree_command(args, alignment, prefix):
    cmd = [
        str(args.iqtree),
        "-s",
        str(alignment),
        "-m",
        args.model,
        "--satute",
        "--satute-alpha",
        str(args.alpha),
        "--prefix",
        str(prefix),
        "-T",
        str(args.threads),
        "--redo",
    ]
    if args.seed is not None:
        cmd.extend(["--seed", str(args.seed)])
    if args.start_tree:
        cmd.extend(["-te", str(args.start_tree)])
        if args.fix_branch_lengths:
            cmd.append("-blfix")
    if args.support_mode in {"alrt", "both"}:
        cmd.extend(["--alrt", str(args.support_replicates)])
    if args.support_mode in {"ufboot", "both"}:
        cmd.extend(["-B", str(args.support_replicates)])
        if args.boot_trees:
            cmd.append("--boot-trees")
    if args.iqtree_arg:
        for extra in args.iqtree_arg:
            cmd.extend(extra.split())
    if args.quiet:
        cmd.append("--quiet")
    return cmd


def write_detail_rows(dataset_name, alignment, prefix, detail_writer, support_mode):
    report = parse_iqtree_report(Path(str(prefix) + ".iqtree"))
    sat_rows = parse_sat_stat(Path(str(prefix) + ".sat.stat"))
    tree_path = Path(str(prefix) + ".treefile")
    if not tree_path.exists():
        tree_path = Path(str(prefix) + ".sat.tree")
    support_by_split, tree_branch_lengths, taxa = tree_support_by_split(tree_path, support_mode)

    for row in sat_rows:
        split = split_from_sat_row(row, taxa)
        support = support_by_split.get(split, {}) if split is not None else {}
        tree_branch_length = tree_branch_lengths.get(split, "") if split is not None else ""
        detail_writer.writerow(
            {
                "dataset": dataset_name,
                "alignment": str(alignment),
                "prefix": str(prefix),
                "taxa": len(taxa),
                "log_likelihood": report["log_likelihood"],
                "unconstrained_log_likelihood": report["unconstrained_log_likelihood"],
                "aic": report["aic"],
                "aicc": report["aicc"],
                "bic": report["bic"],
                "branch_id": row.get("ID", ""),
                "formula": row.get("Formula", "dominant"),
                "branch_length": row.get("Length", tree_branch_length),
                "effective_branch_length": row.get("EffectiveLength", ""),
                "rate_category": row.get("RateCategory", "pooled"),
                "rate_multiplier": row.get("RateMultiplier", ""),
                "satC": row.get("satC", ""),
                "satVar": row.get("satVar", ""),
                "satSE": row.get("satSE", ""),
                "satZ": row.get("satZ", ""),
                "satP": row.get("satP", ""),
                "decision": row.get("Decision", ""),
                "decision_taxon_bonf": row.get("DecisionTaxonBonf", ""),
                "alpha": row.get("Alpha", ""),
                "alpha_adjust": row.get("AlphaTaxonBonf", ""),
                "left_taxa": row.get("LeftTaxa", ""),
                "right_taxa": row.get("RightTaxa", ""),
                "split": row.get("Split", ""),
                "support_label": support.get("support_label", ""),
                "sh_alrt": support.get("sh_alrt", ""),
                "ufboot": support.get("ufboot", ""),
                "tree_split_joined": "1" if support else "0",
                "support_present": "1" if support.get("support_label", "") else "0",
            }
        )


def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def summarize_detail(detail_path, summary_path):
    groups = defaultdict(list)
    with open(detail_path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            groups[(row["dataset"], row["formula"])].append(row)

    fieldnames = [
        "dataset",
        "formula",
        "branches",
        "informative",
        "saturated",
        "bonferroni_informative",
        "bonferroni_saturated",
        "median_satZ",
        "median_satP",
        "high_support_branches",
        "high_support_saturated",
        "low_support_branches",
        "low_support_saturated",
        "log_likelihood",
        "aic",
        "bic",
    ]
    with open(summary_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for (dataset, formula), rows in sorted(groups.items()):
            satz = [to_float(row["satZ"]) for row in rows if math.isfinite(to_float(row["satZ"]))]
            satp = [to_float(row["satP"]) for row in rows if math.isfinite(to_float(row["satP"]))]
            high = []
            low = []
            for row in rows:
                sh = to_float(row.get("sh_alrt"))
                uf = to_float(row.get("ufboot"))
                values = [value for value in (sh, uf) if math.isfinite(value)]
                if not values:
                    continue
                if min(values) >= 95.0:
                    high.append(row)
                else:
                    low.append(row)
            writer.writerow(
                {
                    "dataset": dataset,
                    "formula": formula,
                    "branches": len(rows),
                    "informative": sum(1 for row in rows if row["decision"] == "informative"),
                    "saturated": sum(1 for row in rows if row["decision"] == "saturated"),
                    "bonferroni_informative": sum(1 for row in rows if row["decision_taxon_bonf"] == "informative"),
                    "bonferroni_saturated": sum(1 for row in rows if row["decision_taxon_bonf"] == "saturated"),
                    "median_satZ": statistics.median(satz) if satz else "",
                    "median_satP": statistics.median(satp) if satp else "",
                    "high_support_branches": len(high),
                    "high_support_saturated": sum(1 for row in high if row["decision"] == "saturated"),
                    "low_support_branches": len(low),
                    "low_support_saturated": sum(1 for row in low if row["decision"] == "saturated"),
                    "log_likelihood": rows[0]["log_likelihood"] if rows else "",
                    "aic": rows[0]["aic"] if rows else "",
                    "bic": rows[0]["bic"] if rows else "",
                }
            )


def dataset_name_for(path, used):
    stem = Path(path).stem
    name = stem
    counter = 2
    while name in used:
        name = f"{stem}_{counter}"
        counter += 1
    used.add(name)
    return name


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Optimize IQ-TREE trees for one or more alignments, run native SatuTe, "
            "and join likelihood, branch support, and saturation rows by tree split."
        )
    )
    parser.add_argument("alignments", nargs="+", help="Input sequence alignments accepted by IQ-TREE.")
    parser.add_argument("--iqtree", default="./build/iqtree3")
    parser.add_argument("--outdir", default="/tmp/iqtree-satute-dataset-analysis")
    parser.add_argument("--model", default="MFP", help="IQ-TREE model string, e.g. MFP, JC, GTR+G4, LG+G4.")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--support-mode", choices=["none", "alrt", "ufboot", "both"], default="alrt")
    parser.add_argument("--support-replicates", type=int, default=1000)
    parser.add_argument("--boot-trees", action="store_true")
    parser.add_argument("--start-tree", default="", help="Optional fixed starting/reference tree passed with -te.")
    parser.add_argument("--fix-branch-lengths", action="store_true", help="Add -blfix when --start-tree is used.")
    parser.add_argument("--seed", type=int)
    parser.add_argument("--iqtree-arg", action="append", default=[], help="Extra IQ-TREE argument string; may be repeated.")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    run_dir = outdir / "runs"
    run_dir.mkdir(parents=True, exist_ok=True)
    detail_path = outdir / "dataset_tree_satute_detail.tsv"
    summary_path = outdir / "dataset_tree_satute_summary.tsv"
    manifest_path = outdir / "manifest.json"

    fieldnames = [
        "dataset",
        "alignment",
        "prefix",
        "taxa",
        "log_likelihood",
        "unconstrained_log_likelihood",
        "aic",
        "aicc",
        "bic",
        "branch_id",
        "formula",
        "branch_length",
        "effective_branch_length",
        "rate_category",
        "rate_multiplier",
        "satC",
        "satVar",
        "satSE",
        "satZ",
        "satP",
        "decision",
        "decision_taxon_bonf",
        "alpha",
        "alpha_adjust",
        "left_taxa",
        "right_taxa",
        "split",
        "support_label",
        "sh_alrt",
        "ufboot",
        "tree_split_joined",
        "support_present",
    ]

    commands = []
    used_names = set()
    with open(detail_path, "w", encoding="utf-8", newline="") as detail_handle:
        writer = csv.DictWriter(detail_handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for alignment_text in args.alignments:
            alignment = Path(alignment_text)
            if not alignment.exists():
                raise FileNotFoundError(alignment)
            dataset = dataset_name_for(alignment, used_names)
            prefix = run_dir / dataset / dataset
            prefix.parent.mkdir(parents=True, exist_ok=True)
            cmd = build_iqtree_command(args, alignment, prefix)
            commands.append(cmd)
            run_command(cmd, prefix.parent / "iqtree_command.log")
            write_detail_rows(dataset, alignment, prefix, writer, args.support_mode)

    summarize_detail(detail_path, summary_path)
    manifest = {
        "iqtree": str(args.iqtree),
        "model": args.model,
        "alpha": args.alpha,
        "support_mode": args.support_mode,
        "support_replicates": args.support_replicates,
        "commands": commands,
        "detail": str(detail_path),
        "summary": str(summary_path),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {detail_path}")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
