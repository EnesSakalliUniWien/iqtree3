#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
ROOT_DIR="${ROOT_DIR:-$(cd "${PROJECT_DIR}/../.." && pwd)}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
IQTREE_BIN="${IQTREE_BIN:-${ROOT_DIR}/build/iqtree3}"
SEQGEN_BIN="${SEQGEN_BIN:?Set SEQGEN_BIN to the Seq-Gen executable path}"
EVONAPS_BRANCH_LENGTHS="${EVONAPS_BRANCH_LENGTHS:-}"
OUTDIR="${OUTDIR:-/work/satute-head-to-head}"
REPS="${REPS:-1000}"
SHARD_INDEX="${SHARD_INDEX:-${AWS_BATCH_JOB_ARRAY_INDEX:-0}}"
SHARD_COUNT="${SHARD_COUNT:-1}"
SIMULATION_MODELS="${SIMULATION_MODELS:-JC}"
EVALUATION_MODELS="${EVALUATION_MODELS:-${SIMULATION_MODELS}}"
MODEL_PAIRS="${MODEL_PAIRS:-}"
SCENARIO_SET="${SCENARIO_SET:-fig2}"
TREE_CASES="${TREE_CASES:-five_external}"
RESUME="${RESUME:-0}"

if [[ "${TREE_CASES}" == *sixteen_internal* && ! -s "${EVONAPS_BRANCH_LENGTHS}" ]]; then
  echo "sixteen_internal requires EVONAPS_BRANCH_LENGTHS pointing to the empirical EvoNAPS branch-length TSV/CSV" >&2
  exit 1
fi

mkdir -p "${OUTDIR}"

cmd=(
  "${PYTHON_BIN}" "${PROJECT_DIR}/experiments/002_relative_weighting/run.py"
  --iqtree "${IQTREE_BIN}" \
  --seqgen "${SEQGEN_BIN}" \
  --simulator seq-gen \
  --outdir "${OUTDIR}/shard-${SHARD_INDEX}" \
  --reps "${REPS}" \
  --paper-grid \
  --tree-cases "${TREE_CASES}" \
  --simulation-models "${SIMULATION_MODELS}" \
  --evaluation-models "${EVALUATION_MODELS}" \
  --scenario-set "${SCENARIO_SET}" \
  --shard-index "${SHARD_INDEX}" \
  --shard-count "${SHARD_COUNT}"
)

if [[ -n "${EVONAPS_BRANCH_LENGTHS}" ]]; then
  cmd+=(--evonaps-branch-lengths "${EVONAPS_BRANCH_LENGTHS}")
fi

if [[ -n "${MODEL_PAIRS}" ]]; then
  cmd+=(--model-pairs "${MODEL_PAIRS}")
fi

if [[ "${RESUME}" == "1" || "${RESUME}" == "true" ]]; then
  cmd+=(--resume)
fi

"${cmd[@]}"
