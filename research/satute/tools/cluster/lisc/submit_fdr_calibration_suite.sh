#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: submit_fdr_calibration_suite.sh exact|null-mixed [--dry-run] [sbatch options...]" >&2
  exit 2
fi
DESIGN_NAME="$1"
shift
DRY_RUN=0
if [[ "${1:-}" == "--dry-run" ]]; then
  DRY_RUN=1
  shift
fi

ROOT_DIR="${ROOT_DIR:-$(git rev-parse --show-toplevel)}"
PROJECT_DIR="${PROJECT_DIR:-${ROOT_DIR}/research/satute}"
SOURCE_ID="${SOURCE_ID:-$(git -C "${ROOT_DIR}" rev-parse --short=12 HEAD)}"
RUN_DATE="${RUN_DATE:-$(date -u +%Y%m%d)}"
export ROOT_DIR PROJECT_DIR SOURCE_ID
export WORK_DIR="${WORK_DIR:-/lisc/data/scratch/menche/sakalli/iq-tree-satute/full-head-to-head}"
export IQTREE_BIN="${IQTREE_BIN:-${WORK_DIR}/iqtree3-build-${SOURCE_ID}/iqtree3}"
export SHARD_COUNT="${SHARD_COUNT:-1000}"
export MODELS="${MODELS:-JC,GTR_PF06346}"
export SITE_LENGTHS="${SITE_LENGTHS:-100,1000}"
export BRANCH_LENGTHS="${BRANCH_LENGTHS:-4,8}"
export TREE_CASES="${TREE_CASES:-five_external,sixteen_internal}"

case "${DESIGN_NAME}" in
  exact)
    export DESIGN=exact_null REPS="${REPS:-5000}" TRAINING_REPS="${TRAINING_REPS:-2500}"
    export RUN_NAME="${RUN_NAME:-fdr-exact-null-reps${REPS}-${SOURCE_ID}-${RUN_DATE}}"
    ;;
  null-mixed)
    export DESIGN=mixed_null REPS="${REPS:-2000}" TRAINING_REPS=0
    export RUN_NAME="${RUN_NAME:-fdr-mixed-null-reps${REPS}-${SOURCE_ID}-${RUN_DATE}}"
    ;;
  *)
    echo "Unknown calibration design: ${DESIGN_NAME}" >&2
    exit 2
    ;;
esac

printf 'design=%s\nrun_name=%s\nreps=%s\ntraining_reps=%s\nmodels=%s\nsites=%s\nbranches=%s\ntrees=%s\nshards=%s\n' \
  "${DESIGN}" "${RUN_NAME}" "${REPS}" "${TRAINING_REPS}" "${MODELS}" \
  "${SITE_LENGTHS}" "${BRANCH_LENGTHS}" "${TREE_CASES}" "${SHARD_COUNT}"
if [[ "${DRY_RUN}" == "1" ]]; then
  exit 0
fi
sbatch --parsable "$@" "${PROJECT_DIR}/tools/cluster/lisc/submit_fdr_calibration.slurm"
