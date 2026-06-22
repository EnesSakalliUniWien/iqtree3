#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
SHARD_ROOT="${SHARD_ROOT:-/work/satute-head-to-head}"
OUTDIR="${OUTDIR:-${ROOT_DIR}/doc/satute-wiki/results/full-run}"
FIGURE_DIR="${FIGURE_DIR:-${ROOT_DIR}/doc/satute-wiki/figures/full-run}"
EXPECTED_REPS="${EXPECTED_REPS:-1000}"
TREE_CASES="${TREE_CASES:-five_external,sixteen_internal}"
PLOT_MODELS="${PLOT_MODELS:-JC}"

mapfile -t DETAILS < <(find "${SHARD_ROOT}" -path '*/head_to_head_detail.tsv' -type f | sort)
if [[ "${#DETAILS[@]}" -eq 0 ]]; then
  echo "No shard detail TSV files found under ${SHARD_ROOT}" >&2
  exit 1
fi

mkdir -p "${OUTDIR}" "${FIGURE_DIR}"

"${PYTHON_BIN}" "${ROOT_DIR}/doc/satute-wiki/combine_head_to_head_shards.py" \
  --inputs "${DETAILS[@]}" \
  --outdir "${OUTDIR}"

"${PYTHON_BIN}" "${ROOT_DIR}/doc/satute-wiki/verify_head_to_head_calculations.py" \
  --detail "${OUTDIR}/head_to_head_detail.tsv"

IFS=',' read -r -a plot_models <<< "${PLOT_MODELS}"
for model in "${plot_models[@]}"; do
  "${PYTHON_BIN}" "${ROOT_DIR}/doc/satute-wiki/plot_head_to_head_full.py" \
    --summary "${OUTDIR}/head_to_head_summary.tsv" \
    --outdir "${FIGURE_DIR}" \
    --expected-reps "${EXPECTED_REPS}" \
    --tree-cases "${TREE_CASES}" \
    --model "${model}"
done

echo "Combined results: ${OUTDIR}"
echo "Full-run figures: ${FIGURE_DIR}"
