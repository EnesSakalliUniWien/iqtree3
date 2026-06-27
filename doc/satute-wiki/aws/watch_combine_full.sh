#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
SHARD_ROOT="${SHARD_ROOT:-/work/satute-head-to-head}"
OUTDIR="${OUTDIR:-${ROOT_DIR}/doc/satute-wiki/results/full-run}"
FIGURE_DIR="${FIGURE_DIR:-${ROOT_DIR}/doc/satute-wiki/figures/full-run}"
EXPECTED_REPS="${EXPECTED_REPS:-1000}"
TREE_CASES="${TREE_CASES:-five_external}"
PLOT_MODEL_PAIRS="${PLOT_MODEL_PAIRS:-JC:JC}"
PLOT_SCENARIO_SET="${PLOT_SCENARIO_SET:-fig2}"
POLL_SECONDS="${POLL_SECONDS:-60}"

worker_count() {
  ps -eo cmd | grep -F "head_to_head_satute_simulations.py" | grep -v grep | wc -l | tr -d ' '
}

detail_rows() {
  find "${SHARD_ROOT}" -name head_to_head_detail.tsv -print0 2>/dev/null \
    | xargs -0 awk 'FNR>1{n++} END{print n+0}' 2>/dev/null || echo 0
}

while true; do
  workers="$(worker_count)"
  rows="$(detail_rows)"
  printf '%s workers=%s rows=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${workers}" "${rows}"
  if [[ "${workers}" == "0" ]]; then
    break
  fi
  sleep "${POLL_SECONDS}"
done

SHARD_ROOT="${SHARD_ROOT}" \
OUTDIR="${OUTDIR}" \
FIGURE_DIR="${FIGURE_DIR}" \
EXPECTED_REPS="${EXPECTED_REPS}" \
TREE_CASES="${TREE_CASES}" \
PLOT_MODEL_PAIRS="${PLOT_MODEL_PAIRS}" \
PLOT_SCENARIO_SET="${PLOT_SCENARIO_SET}" \
"${ROOT_DIR}/doc/satute-wiki/aws/combine_verify_plot_full.sh"
