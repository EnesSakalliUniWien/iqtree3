#!/bin/bash

set -euo pipefail

IQTREE=${1:-./build/iqtree3}
OUTROOT=${2:-/tmp/iqtree-satute-phase1-checks}
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
GTR_REFERENCE_MODEL='GTR{1,2,1,1,2,1}+F{0.30,0.20,0.20,0.30}'
GTR_GAMMA_REFERENCE_MODEL='GTR{1,2,1,1,2,1}+F{0.30,0.20,0.20,0.30}+G4{0.5}'
JC_INVAR_GAMMA_REFERENCE_MODEL='JC+I{0.1}+G4{0.5}'

if [ ! -x "$IQTREE" ]; then
    echo "Cannot execute IQ-TREE binary: $IQTREE" >&2
    echo "Usage: $0 [path/to/iqtree3] [output-root]" >&2
    exit 2
fi

PYTHON_WITH_NUMPY=${PYTHON_BIN:-$(command -v python3 || true)}

if [ -z "${PYTHON_WITH_NUMPY:-}" ]; then
    echo "Cannot find python3 for formula comparison" >&2
    exit 2
fi

if ! "$PYTHON_WITH_NUMPY" - <<'PY'
import numpy
PY
then
    echo "Python runtime lacks NumPy, required by compare_formulas.py: $PYTHON_WITH_NUMPY" >&2
    exit 2
fi

mkdir -p "$OUTROOT"

echo "Running SatuTe phase-1 checks"
echo "  IQ-TREE: $IQTREE"
echo "  Output:  $OUTROOT"
echo "  Python:  $PYTHON_WITH_NUMPY"

echo
echo "[1/5] Simulated JC/GTR native smoke checks"
"$SCRIPT_DIR/test_satute_simulated.sh" "$IQTREE" "$OUTROOT/simulated-smoke"

echo
echo "[2/5] IQ-TREE-owned rate-category checks with independent Python references"
"$SCRIPT_DIR/test_satute_rate_categories.py" "$IQTREE" "$OUTROOT/rate-categories"

echo
echo "[3/5] Formula comparison against Python reference"
"$PYTHON_WITH_NUMPY" "$SCRIPT_DIR/compare_formulas.py" \
    --iqtree "$IQTREE" \
    --outdir "$OUTROOT/formula-comparison"

echo
echo "[4/5] Reusable Python reference CLI check"
REFERENCE_OUT="$OUTROOT/formula-comparison/gtr_5.00_reference.tsv"
"$PYTHON_WITH_NUMPY" "$SCRIPT_DIR/satute_reference.py" \
    --alignment "$OUTROOT/formula-comparison/gtr_5.00.fa" \
    --tree "$OUTROOT/formula-comparison/gtr_5.00_sat.sat.tree" \
    --model "$GTR_REFERENCE_MODEL" \
    --split C,D \
    --out "$REFERENCE_OUT" \
    --compare-sat-stat "$OUTROOT/formula-comparison/gtr_5.00_sat.sat.stat"

RATE_REFERENCE_OUT="$OUTROOT/rate-categories/gtr_gamma/reference_cli_rate.tsv"
"$PYTHON_WITH_NUMPY" "$SCRIPT_DIR/satute_reference.py" \
    --alignment "$OUTROOT/rate-categories/sim.fa" \
    --tree "$OUTROOT/rate-categories/gtr_gamma/example.sat.tree" \
    --model "$GTR_GAMMA_REFERENCE_MODEL" \
    --rate-file "$OUTROOT/rate-categories/gtr_gamma/example.rate" \
    --iqtree-report "$OUTROOT/rate-categories/gtr_gamma/example.iqtree" \
    --out "$RATE_REFERENCE_OUT" \
    --compare-sat-stat "$OUTROOT/rate-categories/gtr_gamma/example.sat.stat"

INVAR_RATE_REFERENCE_OUT="$OUTROOT/rate-categories/invar_gamma/reference_cli_rate.tsv"
"$PYTHON_WITH_NUMPY" "$SCRIPT_DIR/satute_reference.py" \
    --alignment "$OUTROOT/rate-categories/sim.fa" \
    --tree "$OUTROOT/rate-categories/invar_gamma/example.sat.tree" \
    --model "$JC_INVAR_GAMMA_REFERENCE_MODEL" \
    --rate-file "$OUTROOT/rate-categories/invar_gamma/example.rate" \
    --iqtree-report "$OUTROOT/rate-categories/invar_gamma/example.iqtree" \
    --out "$INVAR_RATE_REFERENCE_OUT" \
    --compare-sat-stat "$OUTROOT/rate-categories/invar_gamma/example.sat.stat"

"$PYTHON_WITH_NUMPY" "$SCRIPT_DIR/../exact/test_weighted_likelihood.py"

BAD_SPLIT_LOG="$OUTROOT/formula-comparison/bad_split_reference.log"
if "$PYTHON_WITH_NUMPY" "$SCRIPT_DIR/satute_reference.py" \
    --alignment "$OUTROOT/formula-comparison/gtr_5.00.fa" \
    --tree "$OUTROOT/formula-comparison/gtr_5.00_sat.sat.tree" \
    --model "$GTR_REFERENCE_MODEL" \
    --split A,C \
    --compare-sat-stat "$OUTROOT/formula-comparison/gtr_5.00_sat.sat.stat" \
    >"$BAD_SPLIT_LOG" 2>&1
then
    echo "Reusable Python reference CLI accepted invalid split A,C" >&2
    exit 1
fi
if ! grep -q "not a branch in the tree" "$BAD_SPLIT_LOG"; then
    echo "Reusable Python reference CLI rejected invalid split without expected message" >&2
    cat "$BAD_SPLIT_LOG" >&2
    exit 1
fi
echo "Reusable Python reference CLI rejected invalid split A,C"

echo
echo "[5/5] Generalized larger-tree Python reference check"
"$PYTHON_WITH_NUMPY" "$SCRIPT_DIR/test_satute_reference_general.py" \
    "$IQTREE" \
    "$OUTROOT/general-reference"

echo
echo "SatuTe phase-1 checks passed; outputs are in $OUTROOT"
