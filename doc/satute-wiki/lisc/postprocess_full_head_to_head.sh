#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
WORK_DIR="${WORK_DIR:-/lisc/data/scratch/menche/sakalli/iq-tree-satute/full-head-to-head}"
RUN_NAME="${RUN_NAME:-paper-fig2-jc-reps1000}"
OUT_ROOT="${OUT_ROOT:-${WORK_DIR}/${RUN_NAME}}"
LISC_PYTHON_MODULES="${LISC_PYTHON_MODULES:-SciPy-bundle/2025.07-gfbf-2025b}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
EXPECTED_REPS="${EXPECTED_REPS:-1000}"
TREE_CASES="${TREE_CASES:-five_external,sixteen_internal}"
PLOT_MODEL_PAIRS="${PLOT_MODEL_PAIRS:-JC:JC}"
PLOT_SCENARIO_SET="${PLOT_SCENARIO_SET:-fig2}"

if command -v module >/dev/null 2>&1 && [[ -n "${LISC_PYTHON_MODULES}" ]]; then
  set +u
  # shellcheck disable=SC2086
  module load ${LISC_PYTHON_MODULES}
  set -u
fi

SHARD_ROOT="${OUT_ROOT}/shards" \
OUTDIR="${OUT_ROOT}/merged" \
FIGURE_DIR="${OUT_ROOT}/figures" \
EXPECTED_REPS="${EXPECTED_REPS}" \
TREE_CASES="${TREE_CASES}" \
PLOT_MODEL_PAIRS="${PLOT_MODEL_PAIRS}" \
PLOT_SCENARIO_SET="${PLOT_SCENARIO_SET}" \
PYTHON_BIN="${PYTHON_BIN}" \
ROOT_DIR="${ROOT_DIR}" \
  "${ROOT_DIR}/doc/satute-wiki/aws/combine_verify_plot_full.sh"

echo "Run root: ${OUT_ROOT}"
echo "Run archives: ${OUT_ROOT}/run-archives"
echo "Merged TSVs: ${OUT_ROOT}/merged"
echo "Figures: ${OUT_ROOT}/figures"
