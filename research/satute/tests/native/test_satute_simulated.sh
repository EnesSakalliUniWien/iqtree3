#!/bin/bash

set -euo pipefail

IQTREE=${1:-build/iqtree3}
OUTDIR=${2:-/tmp/iqtree-satute-sim-smoke}
PYTHON_BIN=${PYTHON_BIN:-$(command -v python3 || true)}

if [ ! -x "$IQTREE" ]; then
    echo "Cannot execute IQ-TREE binary: $IQTREE" >&2
    echo "Usage: $0 [path/to/iqtree3] [output-dir]" >&2
    exit 2
fi
if [ -z "${PYTHON_BIN:-}" ]; then
    echo "Cannot find python3 for SatuTe FDR verification" >&2
    exit 2
fi

rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

run_case() {
    local label=$1
    local model=$2
    local branch_length=$3
    local expected=$4
    local seed=$5
    local tree="$OUTDIR/${label}.tree"
    local sim_prefix="$OUTDIR/${label}"
    local sat_prefix="$OUTDIR/${label}_sat"
    local stat_file="${sat_prefix}.sat.stat"

    printf '((A:0.05,B:0.05):%s,C:0.05,D:0.05);\n' "$branch_length" > "$tree"

    "$IQTREE" --alisim "$sim_prefix" -t "$tree" -m "$model" --length 5000 \
        --seed "$seed" -af fasta --quiet
    "$IQTREE" -s "${sim_prefix}.fa" -te "$tree" -m "$model" --satute \
        --prefix "$sat_prefix" -T 1 --redo --quiet

    local row decision pvalue zscore
    row=$(awk 'BEGIN{FS=OFS="\t"}
        $1 == "ID" {
            for (i = 1; i <= NF; i++) h[$i] = i
            next
        }
        $1 ~ /^[0-9]+$/ && $h["Formula"] == "dominant" && $h["RateCategory"] == "pooled" && $h["Split"] == "A,B" {
            print
            exit
        }' "$stat_file")
    decision=$(printf '%s\n' "$row" | awk 'BEGIN{FS="\t"} {print $17}')
    zscore=$(printf '%s\n' "$row" | awk 'BEGIN{FS="\t"} {print $13}')
    pvalue=$(printf '%s\n' "$row" | awk 'BEGIN{FS="\t"} {print $14}')

    printf '%s\tmodel=%s\tbranch=%s\tz=%s\tp=%s\tdecision=%s\n' \
        "$label" "$model" "$branch_length" "$zscore" "$pvalue" "$decision"

    if [ "$decision" != "$expected" ]; then
        echo "Expected $expected for $label, got $decision" >&2
        echo "$row" >&2
        exit 1
    fi
}

verify_fdr() {
    local stat_file=$1
    local alpha=$2
    "$PYTHON_BIN" - "$stat_file" "$alpha" <<'PY'
import csv
import math
import sys
from collections import defaultdict

path = sys.argv[1]
alpha = float(sys.argv[2])
rows = []
with open(path, encoding="utf-8") as handle:
    reader = csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t")
    rows = list(reader)

families = defaultdict(list)
for row_index, row in enumerate(rows):
    formula = row["Formula"]
    if formula not in {"dominant", "eigenvalue_weighted"}:
        raise AssertionError(f"unexpected formula in {path}: {formula}")
    if row["RateCategory"] != "pooled":
        if row["FDR_BY"] != "NA" or row["DecisionFDR"] != "not_tested":
            raise AssertionError(f"category row entered FDR family in {path}: {row}")
        continue
    families[formula].append((float(row["satP"]), row_index, row))

if set(families) != {"dominant", "eigenvalue_weighted"}:
    raise AssertionError(f"missing separate FDR family in {path}: {sorted(families)}")
if len(families["dominant"]) != len(families["eigenvalue_weighted"]):
    raise AssertionError(f"formula FDR family sizes differ in {path}")

for formula, entries in families.items():
    entries.sort(key=lambda item: (item[0], item[1]))
    family_size = len(entries)
    harmonic = sum(1.0 / rank for rank in range(1, family_size + 1))
    expected = [None] * family_size
    running = 1.0
    for offset in range(family_size - 1, -1, -1):
        rank = offset + 1
        running = min(running, min(1.0, entries[offset][0] * family_size * harmonic / rank))
        expected[offset] = running
    for entry, expected_fdr in zip(entries, expected):
        observed = float(entry[2]["FDR_BY"])
        if not math.isclose(observed, expected_fdr, rel_tol=0.0, abs_tol=2e-9):
            raise AssertionError(
                f"{path} {formula}: observed BY {observed} != expected {expected_fdr}"
            )
        expected_decision = "informative" if expected_fdr <= alpha else "saturated"
        if entry[2]["DecisionFDR"] != expected_decision:
            raise AssertionError(
                f"{path} {formula}: decision {entry[2]['DecisionFDR']} != {expected_decision}"
            )

print(
    f"fdr_verified\tfile={path}\tdominant={len(families['dominant'])}"
    f"\teigenvalue_weighted={len(families['eigenvalue_weighted'])}"
)
PY
}

GTR_MODEL='GTR{1,2,1,1,2,1}+F{0.30,0.20,0.20,0.30}'

echo -e "case\tmodel\tbranch\tz\tp\tdecision"
run_case jc_informative JC 0.50 informative 11
run_case jc_saturated JC 8.00 saturated 800
run_case gtr_informative "$GTR_MODEL" 0.50 informative 50
run_case gtr_saturated "$GTR_MODEL" 8.00 saturated 800

verify_fdr "$OUTDIR/jc_informative_sat.sat.stat" 0.05
verify_fdr "$OUTDIR/jc_saturated_sat.sat.stat" 0.05
verify_fdr "$OUTDIR/gtr_informative_sat.sat.stat" 0.05
verify_fdr "$OUTDIR/gtr_saturated_sat.sat.stat" 0.05

edge_file="$OUTDIR/jc_informative_edges.txt"
subset_prefix="$OUTDIR/jc_informative_subset"
selected_edge=$(awk 'BEGIN{FS="\t"}
    $1 == "ID" {
        for (i = 1; i <= NF; i++) h[$i] = i
        next
    }
    $1 ~ /^[0-9]+$/ && $h["Formula"] == "dominant" && $h["RateCategory"] == "pooled" && $h["Split"] == "A,B" {
        print $h["ID"]
        exit
    }' "$OUTDIR/jc_informative_sat.sat.stat")

if [ -z "$selected_edge" ]; then
    echo "Could not find A,B branch ID in jc_informative SatuTe output" >&2
    exit 1
fi

printf '%s\n' "$selected_edge" > "$edge_file"
"$IQTREE" -s "$OUTDIR/jc_informative.fa" -te "$OUTDIR/jc_informative.tree" -m JC \
    --satute --satute-edges "$edge_file" --satute-alpha 0.01 \
    --prefix "$subset_prefix" -T 1 --redo --quiet

subset_summary=$(awk -v edge="$selected_edge" 'BEGIN{FS="\t"; rows=0; bad=0; alpha_bad=0; dominant=0; weighted=0; unexpected=0}
    $1 == "ID" {
        for (i = 1; i <= NF; i++) h[$i] = i
        next
    }
    $1 ~ /^[0-9]+$/ {
        rows++
        if ($h["ID"] != edge) bad++
        if ($h["Alpha"] != "0.01") alpha_bad++
        if ($h["Formula"] == "dominant") dominant++
        else if ($h["Formula"] == "eigenvalue_weighted") weighted++
        else unexpected++
    }
    END {
        printf "%d\t%d\t%d\t%d\t%d\t%d", rows, bad, alpha_bad, dominant, weighted, unexpected
    }' "$subset_prefix.sat.stat")
subset_rows=$(printf '%s\n' "$subset_summary" | awk 'BEGIN{FS="\t"} {print $1}')
subset_bad=$(printf '%s\n' "$subset_summary" | awk 'BEGIN{FS="\t"} {print $2}')
subset_alpha_bad=$(printf '%s\n' "$subset_summary" | awk 'BEGIN{FS="\t"} {print $3}')
subset_dominant=$(printf '%s\n' "$subset_summary" | awk 'BEGIN{FS="\t"} {print $4}')
subset_weighted=$(printf '%s\n' "$subset_summary" | awk 'BEGIN{FS="\t"} {print $5}')
subset_unexpected=$(printf '%s\n' "$subset_summary" | awk 'BEGIN{FS="\t"} {print $6}')

printf 'jc_subset\tmodel=JC\tedge=%s\trows=%s\talpha=0.01\n' "$selected_edge" "$subset_rows"

if [ "$subset_rows" != "2" ] || [ "$subset_bad" != "0" ] || [ "$subset_alpha_bad" != "0" ] || \
    [ "$subset_dominant" != "1" ] || [ "$subset_weighted" != "1" ] || \
    [ "$subset_unexpected" != "0" ]; then
    echo "SatuTe edge-subset check failed: rows=$subset_rows bad_id=$subset_bad bad_alpha=$subset_alpha_bad dominant=$subset_dominant weighted=$subset_weighted unexpected=$subset_unexpected" >&2
    exit 1
fi

verify_fdr "$subset_prefix.sat.stat" 0.01

if ! awk 'BEGIN{FS="\t"; seen=0; bad=0}
    $1 == "ID" {
        for (i = 1; i <= NF; i++) h[$i] = i
        if (!("InformationFraction" in h) || !("SaturationIndex" in h)) bad++
        next
    }
    $1 ~ /^[0-9]+$/ {
        info = $h["InformationFraction"] + 0
        sat = $h["SaturationIndex"] + 0
        if (info < -1e-12 || info > 1 + 1e-12 || sat < -1e-12 || sat > 1 + 1e-12) bad++
        if ((info + sat - 1) > 1e-9 || (1 - info - sat) > 1e-9) bad++
        if (!seen) { first_info=info; first_sat=sat; seen=1 }
        else if ((info-first_info > 1e-9 || first_info-info > 1e-9) ||
                 (sat-first_sat > 1e-9 || first_sat-sat > 1e-9)) bad++
    }
    END { exit (seen && bad == 0) ? 0 : 1 }
' "$subset_prefix.sat.stat"; then
    echo "SatuTe saturation-scale smoke check failed" >&2
    exit 1
fi

if ! awk 'BEGIN{FS="\t"; checked=0}
    $1 == "ID" { for (i = 1; i <= NF; i++) h[$i] = i; next }
    $1 ~ /^[0-9]+$/ && $h["Formula"] == "dominant" && $h["RateCategory"] == "pooled" && $h["Split"] == "A,B" {
        expected = 1 - exp((-8.0/3.0) * $h["Length"])
        observed = $h["SaturationIndex"] + 0
        difference = observed - expected
        if (difference < 0) difference = -difference
        checked = 1
        exit difference > 1e-8
    }
    END { if (!checked) exit 1 }
' "$subset_prefix.sat.stat"; then
    echo "SatuTe JC closed-form saturation-index check failed" >&2
    exit 1
fi

if ! grep -q 'satIndex=' "$subset_prefix.sat.tree.nex" || ! grep -q 'satInfo=' "$subset_prefix.sat.tree.nex"; then
    echo "SatuTe saturation-scale tree annotations are missing" >&2
    exit 1
fi

if ! grep -q 'satFormula="eigenvalue_weighted"' "$subset_prefix.sat.tree.nex"; then
    echo "SatuTe tree annotations do not identify eigenvalue_weighted as their source formula" >&2
    exit 1
fi
if ! grep -q 'satDominantFDR=' "$subset_prefix.sat.tree.nex" || \
   ! grep -q 'satFDR=' "$subset_prefix.sat.tree.nex" || \
   ! grep -q 'satFDRDecision=' "$subset_prefix.sat.tree.nex"; then
    echo "SatuTe formula-specific FDR tree annotations are missing" >&2
    exit 1
fi

bad_edge_file="$OUTDIR/jc_informative_bad_edges.txt"
bad_edge_prefix="$OUTDIR/jc_informative_bad_subset"
bad_edge_log="$OUTDIR/jc_informative_bad_subset.log"
printf '999999\n' > "$bad_edge_file"
if "$IQTREE" -s "$OUTDIR/jc_informative.fa" -te "$OUTDIR/jc_informative.tree" -m JC \
    --satute --satute-edges "$bad_edge_file" \
    --prefix "$bad_edge_prefix" -T 1 --redo --quiet >"$bad_edge_log" 2>&1
then
    echo "SatuTe edge-subset check accepted an absent branch ID" >&2
    exit 1
fi
if ! grep -q "Requested SatuTe branch IDs not found: 999999" "$bad_edge_log"; then
    echo "SatuTe edge-subset check rejected absent branch ID with unexpected message" >&2
    cat "$bad_edge_log" >&2
    exit 1
fi
printf 'jc_subset_missing\tmodel=JC\tedge=999999\texpected_failure\n'

rooted_tree="$OUTDIR/jc_informative_rooted.tree"
rooted_prefix="$OUTDIR/jc_informative_rooted"
rooted_log="$OUTDIR/jc_informative_rooted.log"
printf '((A:0.05,B:0.05):0.25,(C:0.05,D:0.05):0.25);\n' > "$rooted_tree"
if "$IQTREE" -s "$OUTDIR/jc_informative.fa" -te "$rooted_tree" -m JC \
    --satute --prefix "$rooted_prefix" -T 1 --redo --quiet >"$rooted_log" 2>&1
then
    echo "SatuTe accepted a rooted tree" >&2
    exit 1
fi
if ! grep -q "SatuTe currently supports unrooted trees only" "$rooted_log"; then
    echo "SatuTe rejected a rooted tree with an unexpected message" >&2
    cat "$rooted_log" >&2
    exit 1
fi
printf 'jc_rooted\tmodel=JC\texpected_failure\n'

echo "SatuTe simulated smoke test passed; outputs are in $OUTDIR"
