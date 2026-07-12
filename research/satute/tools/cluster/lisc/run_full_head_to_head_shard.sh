#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
ROOT_DIR="${ROOT_DIR:-$(cd "${PROJECT_DIR}/../.." && pwd)}"
WORK_DIR="${WORK_DIR:-/lisc/data/scratch/menche/sakalli/iq-tree-satute/full-head-to-head}"
RUN_NAME="${RUN_NAME:-paper-fig2-jc-reps1000}"
OUT_ROOT="${OUT_ROOT:-${WORK_DIR}/${RUN_NAME}}"
IQTREE_BIN="${IQTREE_BIN:-${WORK_DIR}/iqtree3-build/iqtree3}"
SEQGEN_BIN="${SEQGEN_BIN:-${WORK_DIR}/seq-gen-1.3.4-source/seq-gen}"
EVONAPS_BRANCH_LENGTHS="${EVONAPS_BRANCH_LENGTHS:-${PROJECT_DIR}/references/evonaps/evonaps_16taxon_branch_lengths.tsv}"
LISC_PYTHON_MODULES="${LISC_PYTHON_MODULES:-SciPy-bundle/2025.07-gfbf-2025b}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
REPS="${REPS:-1000}"
SIMULATION_MODELS="${SIMULATION_MODELS:-JC}"
EVALUATION_MODELS="${EVALUATION_MODELS:-${SIMULATION_MODELS}}"
MODEL_PAIRS="${MODEL_PAIRS:-}"
SCENARIO_SET="${SCENARIO_SET:-fig2}"
TREE_CASES="${TREE_CASES:-five_external,sixteen_internal}"
SHARD_INDEX="${SHARD_INDEX:-${SLURM_ARRAY_TASK_ID:-0}}"
SHARD_COUNT="${SHARD_COUNT:-768}"
ARCHIVE_RUNS="${ARCHIVE_RUNS:-1}"

if command -v module >/dev/null 2>&1 && [[ -n "${LISC_PYTHON_MODULES}" ]]; then
  set +u
  # shellcheck disable=SC2086
  module load ${LISC_PYTHON_MODULES}
  set -u
fi

if [[ ! -x "${IQTREE_BIN}" ]]; then
  echo "Missing IQ-TREE binary: ${IQTREE_BIN}" >&2
  exit 2
fi

if [[ ! -x "${SEQGEN_BIN}" ]]; then
  echo "Missing Seq-Gen binary: ${SEQGEN_BIN}" >&2
  exit 2
fi

if [[ ! -s "${EVONAPS_BRANCH_LENGTHS}" ]]; then
  echo "Missing EvoNAPS branch-length table: ${EVONAPS_BRANCH_LENGTHS}" >&2
  exit 2
fi

SHARD_DIR="${OUT_ROOT}/shards/shard-${SHARD_INDEX}"
mkdir -p "${OUT_ROOT}/logs" "${OUT_ROOT}/run-archives" "${SHARD_DIR}"

model_pair_args=()
if [[ -n "${MODEL_PAIRS}" ]]; then
  model_pair_args=(--model-pairs "${MODEL_PAIRS}")
fi

PYTHONPATH="${PROJECT_DIR}/src${PYTHONPATH:+:${PYTHONPATH}}" "${PYTHON_BIN}" \
  "${PROJECT_DIR}/experiments/002_relative_weighting/run.py" \
  --iqtree "${IQTREE_BIN}" \
  --seqgen "${SEQGEN_BIN}" \
  --simulator seq-gen \
  --evonaps-branch-lengths "${EVONAPS_BRANCH_LENGTHS}" \
  --outdir "${SHARD_DIR}" \
  --reps "${REPS}" \
  --paper-grid \
  --tree-cases "${TREE_CASES}" \
  --simulation-models "${SIMULATION_MODELS}" \
  --evaluation-models "${EVALUATION_MODELS}" \
  "${model_pair_args[@]}" \
  --scenario-set "${SCENARIO_SET}" \
  --shard-index "${SHARD_INDEX}" \
  --shard-count "${SHARD_COUNT}" \
  --resume

if [[ "${ARCHIVE_RUNS}" == "1" || "${ARCHIVE_RUNS}" == "true" ]]; then
  RUNS_DIR="${SHARD_DIR}/runs"
  if [[ -d "${RUNS_DIR}" ]]; then
    ARCHIVE="${OUT_ROOT}/run-archives/runs_shard_${SHARD_INDEX}_job_${SLURM_JOB_ID:-manual}.tar.gz"
    TMP_ARCHIVE="${ARCHIVE}.tmp"
    tar -C "${SHARD_DIR}" -czf "${TMP_ARCHIVE}" runs
    mv "${TMP_ARCHIVE}" "${ARCHIVE}"
    rm -rf "${RUNS_DIR}"
    echo "Archived shard runs: ${ARCHIVE}"
  fi
fi
