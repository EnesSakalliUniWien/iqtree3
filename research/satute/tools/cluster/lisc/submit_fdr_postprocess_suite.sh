#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: submit_fdr_postprocess_suite.sh exact|null-mixed ARRAY_JOB_ID" >&2
  exit 2
fi
DESIGN_NAME="$1"
ARRAY_JOB_ID="$2"

PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
ROOT_DIR="${ROOT_DIR:-$(cd "${PROJECT_DIR}/../.." && pwd)}"
WORK_DIR="${WORK_DIR:-/lisc/data/scratch/menche/sakalli/iq-tree-satute/full-head-to-head}"
SOURCE_ID="${SOURCE_ID:-$(cat "${ROOT_DIR}/SOURCE_COMMIT" 2>/dev/null || git -C "${ROOT_DIR}" rev-parse --short=12 HEAD)}"
SIMULATION_SOURCE_ID="${SIMULATION_SOURCE_ID:-a932abf97b2}"
RUN_DATE="${RUN_DATE:-20260715}"
EXPECTED_SHARDS="${EXPECTED_SHARDS:-1000}"

case "${DESIGN_NAME}" in
  exact)
    EXPECTED_REPS="${EXPECTED_REPS:-5000}"
    TRAINING_REPS="${TRAINING_REPS:-2500}"
    RUN_NAME="${RUN_NAME:-fdr-exact-null-reps5000-${SIMULATION_SOURCE_ID}-${RUN_DATE}}"
    ;;
  null-mixed)
    EXPECTED_REPS="${EXPECTED_REPS:-2000}"
    TRAINING_REPS=0
    RUN_NAME="${RUN_NAME:-fdr-mixed-null-reps2000-${SIMULATION_SOURCE_ID}-${RUN_DATE}}"
    ;;
  *)
    echo "Unknown calibration design: ${DESIGN_NAME}" >&2
    exit 2
    ;;
esac

RUN_ROOT="${RUN_ROOT:-${WORK_DIR}/${RUN_NAME}}"
OUTPUT_ROOT="${OUTPUT_ROOT:-${WORK_DIR}/${RUN_NAME}-postprocessed-${SOURCE_ID}}"
INPUT_LOG_DIR="${INPUT_LOG_DIR:-/lisc/home/user/sakalli/projects/iq-tree-satute-${SIMULATION_SOURCE_ID}/research/satute/artifacts/cluster/lisc/logs}"
export ROOT_DIR PROJECT_DIR WORK_DIR SOURCE_ID ARRAY_JOB_ID EXPECTED_SHARDS
export EXPECTED_REPS TRAINING_REPS RUN_ROOT OUTPUT_ROOT INPUT_LOG_DIR
mkdir -p "${PROJECT_DIR}/artifacts/cluster/lisc/logs"
cd "${ROOT_DIR}"
printf 'design=%s\narray_job=%s\nrun_root=%s\noutput_root=%s\n' \
  "${DESIGN_NAME}" "${ARRAY_JOB_ID}" "${RUN_ROOT}" "${OUTPUT_ROOT}"
if [[ "${DRY_RUN:-0}" == "1" ]]; then
  printf 'source_id=%s\nsimulation_source_id=%s\ninput_log_dir=%s\n' \
    "${SOURCE_ID}" "${SIMULATION_SOURCE_ID}" "${INPUT_LOG_DIR}"
  printf 'postprocess_script=%s\n' \
    "${PROJECT_DIR}/tools/cluster/lisc/postprocess_fdr_calibration.slurm"
  exit 0
fi
sbatch --parsable "${PROJECT_DIR}/tools/cluster/lisc/postprocess_fdr_calibration.slurm"
