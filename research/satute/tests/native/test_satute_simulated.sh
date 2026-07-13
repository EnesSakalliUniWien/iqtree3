#!/bin/bash

set -euo pipefail

IQTREE=${1:-build/iqtree3}
OUTDIR=${2:-/tmp/iqtree-satute-sim-smoke}

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

GTR_MODEL='GTR{1,2,1,1,2,1}+F{0.30,0.20,0.20,0.30}'

echo -e "case\tmodel\tbranch\tz\tp\tdecision"
run_case jc_informative JC 0.50 informative 11
run_case jc_saturated JC 8.00 saturated 800
run_case gtr_informative "$GTR_MODEL" 0.50 informative 50
run_case gtr_saturated "$GTR_MODEL" 8.00 saturated 800

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

subset_summary=$(awk -v edge="$selected_edge" 'BEGIN{FS="\t"; rows=0; bad=0; alpha_bad=0; dominant=0; weighted=0; mixture=0; unexpected=0}
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
        else if ($h["Formula"] == "mixture_likelihood_weighted") mixture++
        else unexpected++
    }
    END {
        printf "%d\t%d\t%d\t%d\t%d\t%d\t%d", rows, bad, alpha_bad, dominant, weighted, mixture, unexpected
    }' "$subset_prefix.sat.stat")
subset_rows=$(printf '%s\n' "$subset_summary" | awk 'BEGIN{FS="\t"} {print $1}')
subset_bad=$(printf '%s\n' "$subset_summary" | awk 'BEGIN{FS="\t"} {print $2}')
subset_alpha_bad=$(printf '%s\n' "$subset_summary" | awk 'BEGIN{FS="\t"} {print $3}')
subset_dominant=$(printf '%s\n' "$subset_summary" | awk 'BEGIN{FS="\t"} {print $4}')
subset_weighted=$(printf '%s\n' "$subset_summary" | awk 'BEGIN{FS="\t"} {print $5}')
subset_mixture=$(printf '%s\n' "$subset_summary" | awk 'BEGIN{FS="\t"} {print $6}')
subset_unexpected=$(printf '%s\n' "$subset_summary" | awk 'BEGIN{FS="\t"} {print $7}')

printf 'jc_subset\tmodel=JC\tedge=%s\trows=%s\talpha=0.01\n' "$selected_edge" "$subset_rows"

if [ "$subset_rows" != "3" ] || [ "$subset_bad" != "0" ] || [ "$subset_alpha_bad" != "0" ] || \
    [ "$subset_dominant" != "1" ] || [ "$subset_weighted" != "1" ] || [ "$subset_mixture" != "1" ] || \
    [ "$subset_unexpected" != "0" ]; then
    echo "SatuTe edge-subset check failed: rows=$subset_rows bad_id=$subset_bad bad_alpha=$subset_alpha_bad dominant=$subset_dominant weighted=$subset_weighted mixture=$subset_mixture unexpected=$subset_unexpected" >&2
    exit 1
fi

if ! awk 'BEGIN{FS="\t"}
    $1 == "ID" { for (i = 1; i <= NF; i++) h[$i] = i; next }
    $1 ~ /^[0-9]+$/ && $h["RateCategory"] == "pooled" {
        if ($h["Formula"] == "eigenvalue_weighted") weighted = $h["satC"]
        if ($h["Formula"] == "mixture_likelihood_weighted") mixture = $h["satC"]
    }
    END {
        if (weighted == "" || mixture == "") exit 1
        difference = (weighted + 0) - (mixture + 0)
        if (difference < 0) difference = -difference
        exit difference > 1e-8
    }
' "$subset_prefix.sat.stat"; then
    echo "Homogeneous JC eigenvalue-weighted and mixture coherence values differ" >&2
    exit 1
fi

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

echo "SatuTe simulated smoke test passed; outputs are in $OUTDIR"
