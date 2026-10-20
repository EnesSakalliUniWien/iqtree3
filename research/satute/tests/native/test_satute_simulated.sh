#!/bin/bash

set -euo pipefail

IQTREE=${1:-build/iqtree3}
OUTDIR=${2:-/tmp/iqtree-satute-sim-smoke}
: "${PYTHON_BIN:?Invoke this adapter through ./satute check all}"

if [ ! -x "$IQTREE" ]; then
    echo "Cannot execute IQ-TREE binary: $IQTREE" >&2
    echo "Usage: $0 [path/to/iqtree3] [output-dir]" >&2
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
# A reversible model unroots the tree, and the root carries no information, so
# the rows must equal those of the unrooted run.
if ! "$IQTREE" -s "$OUTDIR/jc_informative.fa" -te "$rooted_tree" -m JC \
    --satute --prefix "$rooted_prefix" -T 1 --redo --quiet >"$rooted_log" 2>&1
then
    echo "SatuTe rejected a rooted tree under a reversible model" >&2
    cat "$rooted_log" >&2
    exit 1
fi
stat_rows() {
    awk 'BEGIN{FS="\t"} $1 == "ID" { for (i = 1; i <= NF; i++) h[$i] = i; next }
        $1 ~ /^[0-9]+$/ { print $h["Split"], $h["Formula"], $h["RateCategory"], $h["satZ"], $h["SaturationIndex"] }' "$1" | sort
}
if [ "$(stat_rows "$rooted_prefix.sat.stat")" != "$(stat_rows "$OUTDIR/jc_informative_sat.sat.stat")" ]; then
    echo "SatuTe rooted-input JC rows differ from the unrooted run" >&2
    exit 1
fi
printf 'jc_rooted\tmodel=JC\tmatches=unrooted\n'

# Non-reversible models run on IQ-TREE's rooted tree, where the two edges of a
# two-child root are one split: both rows must report the whole root path.
# Equal UNREST rates give JC dynamics, so the index is 1 - exp(-8(t1+t2)/3).
root_tree="$OUTDIR/unrest_root.tree"
root_prefix="$OUTDIR/unrest_root"
root_model='UNREST{1,1,1,1,1,1,1,1,1,1,1}'
printf '((A:0.05,B:0.05):0.5,C:0.05,(D:0.05,E:0.05):0.3);\n' > "$root_tree"
"$IQTREE" --alisim "$root_prefix" -t "$root_tree" -m "$root_model" --length 5000 \
    --seed 11 -af fasta --redo --quiet >/dev/null
"$IQTREE" -s "$root_prefix.fa" -te "$root_tree" -m "$root_model" --satute \
    --prefix "${root_prefix}_sat" -T 1 --redo --quiet >/dev/null

# Exactly one split may span two rows (the root edges), and both must match
# the closed form; a second argument also fixes which split that is.
check_root_path() {
    awk -v want="${2:-}" 'BEGIN{FS="\t"}
        $1 == "ID" { for (i = 1; i <= NF; i++) h[$i] = i; next }
        $1 ~ /^[0-9]+$/ && $h["Formula"] == "dominant" && $h["RateCategory"] == "pooled" {
            split_id = $h["Split"]
            rows[split_id]++
            path_length[split_id] += $h["Length"]
            index_of[split_id, rows[split_id]] = $h["SaturationIndex"] + 0
        }
        END {
            pairs = 0
            for (split_id in rows) {
                if (rows[split_id] != 2) continue
                pairs++
                if (want != "" && split_id != want) exit 1
                expected = 1 - exp((-8.0/3.0) * path_length[split_id])
                for (k = 1; k <= 2; k++) {
                    difference = index_of[split_id, k] - expected
                    if (difference < 0) difference = -difference
                    if (difference > 1e-8) exit 1
                }
            }
            exit (pairs == 1) ? 0 : 1
        }
    ' "$1"
}

if ! check_root_path "${root_prefix}_sat.sat.stat"; then
    echo "SatuTe root-path saturation-index check failed" >&2
    exit 1
fi
printf 'unrest_root_path\tmodel=UNREST\tsplit=root\n'

# A user-rooted tree is accepted for non-reversible models and its root kept:
# rooted on the D,E edge, unlike the root IQ-TREE picks for the unrooted tree.
user_root_tree="$OUTDIR/unrest_user_root.tree"
user_root_prefix="$OUTDIR/unrest_user_root"
printf '(((A:0.05,B:0.05):0.5,C:0.05):0.1,(D:0.05,E:0.05):0.2);\n' > "$user_root_tree"
if ! "$IQTREE" -s "$root_prefix.fa" -te "$user_root_tree" -m "$root_model" --satute \
    --prefix "$user_root_prefix" -T 1 --redo --quiet >"$user_root_prefix.out" 2>&1; then
    echo "SatuTe rejected a user-rooted tree for a non-reversible model" >&2
    cat "$user_root_prefix.out" >&2
    exit 1
fi
if ! check_root_path "$user_root_prefix.sat.stat" "D,E"; then
    echo "SatuTe did not keep the user-supplied root for a non-reversible model" >&2
    exit 1
fi
printf 'unrest_user_root\tmodel=UNREST\tsplit=D,E\n'

# The root-path operator M(t1)^T M(t2) depends on where the root sits on its
# edge; M(t1 + t2) and M(t1) M(t2) do not. With a fixed, strongly
# non-reversible UNREST+G4 and fixed lengths, both root rows must agree in
# every rate category, and moving the root along the same edge must change
# the index (0.5906 at 0.05+0.55 versus 0.5874 at 0.3+0.3).
nonrev_model='UNREST{0.641,0.046,0.028,0.902,0.199,5.114,2.747,0.079,0.038,0.018,0.096}+G4{0.5}'
printf '(((A:0.05,B:0.05):0.5,C:0.05):0.05,(D:0.05,E:0.05):0.55);\n' > "$OUTDIR/nonrev_near.tree"
printf '(((A:0.05,B:0.05):0.5,C:0.05):0.3,(D:0.05,E:0.05):0.3);\n' > "$OUTDIR/nonrev_middle.tree"
for position in near middle; do
    "$IQTREE" -s "$root_prefix.fa" -te "$OUTDIR/nonrev_$position.tree" -m "$nonrev_model" \
        -blfix --satute --prefix "$OUTDIR/nonrev_$position" -T 1 --redo --quiet >/dev/null
done
if ! awk 'BEGIN{FS="\t"}
    FNR == 1 { file++ }
    $1 == "ID" { for (i = 1; i <= NF; i++) h[$i] = i; next }
    $1 ~ /^[0-9]+$/ && $h["Formula"] == "dominant" && $h["Split"] == "D,E" {
        key = file SUBSEP $h["RateCategory"]
        value = $h["SaturationIndex"] + 0
        if (key in first) {
            difference = value - first[key]
            if (difference < 0) difference = -difference
            if (difference > 1e-9) bad++
            paired++
        } else first[key] = value
        if ($h["RateCategory"] == "pooled") pooled[file] = value
    }
    END {
        gap = pooled[1] - pooled[2]
        if (gap < 0) gap = -gap
        exit (bad == 0 && paired == 10 && gap > 1e-4) ? 0 : 1
    }
' "$OUTDIR/nonrev_near.sat.stat" "$OUTDIR/nonrev_middle.sat.stat"; then
    echo "SatuTe non-reversible root-position check failed" >&2
    exit 1
fi
printf 'unrest_g4_root_position\tmodel=UNREST+G4\tsplit=D,E\n'

echo "SatuTe simulated smoke test passed; outputs are in $OUTDIR"
