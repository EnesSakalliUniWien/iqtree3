#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: submit_publication_suite.sh DESIGN [--dry-run] [sbatch options...]

DESIGN is one of:
  gtr       Matched GTR, fixed GTR+G4, and fixed GTR+I+G4
  protein   Matched LG, WAG, JTT, and Q.PFAM with fixed G4 shape
  misspec   GTR data evaluated by GTR, JC, K2P, and F81

The script exports the complete run specification through the process
environment. It deliberately does not use `sbatch --export`, whose comma
separator can truncate model-pair and grid values.
EOF
}

if [[ $# -lt 1 ]]; then
  usage >&2
  exit 2
fi

DESIGN="$1"
shift
DRY_RUN=0
if [[ "${1:-}" == "--dry-run" ]]; then
  DRY_RUN=1
  shift
fi

ROOT_DIR="${ROOT_DIR:-$(git rev-parse --show-toplevel)}"
PROJECT_DIR="${PROJECT_DIR:-${ROOT_DIR}/research/satute}"
SUBMIT_SCRIPT="${PROJECT_DIR}/tools/cluster/lisc/submit_full_head_to_head.slurm"
if [[ ! -s "${SUBMIT_SCRIPT}" ]]; then
  echo "Missing Slurm submission script: ${SUBMIT_SCRIPT}" >&2
  exit 2
fi

SOURCE_ID="${SOURCE_ID:-$(git -C "${ROOT_DIR}" rev-parse --short=12 HEAD)}"
RUN_DATE="${RUN_DATE:-$(date -u +%Y%m%d)}"
export ROOT_DIR PROJECT_DIR
export WORK_DIR="${WORK_DIR:-/lisc/data/scratch/menche/sakalli/iq-tree-satute/full-head-to-head}"
export IQTREE_BIN="${IQTREE_BIN:-${WORK_DIR}/iqtree3-build-${SOURCE_ID}/iqtree3}"
export REPS="${REPS:-1000}"
export SHARD_COUNT="${SHARD_COUNT:-1000}"
export TREE_CASES="${TREE_CASES:-five_external,sixteen_internal}"
export SITE_LENGTHS="${SITE_LENGTHS:-100,250,1000}"
export BRANCH_LENGTHS="${BRANCH_LENGTHS:-4,5,6,7,8,9,10,12}"
export PAPER_GRID=0
export SIMULATOR=alisim
export RESUME_SHARDS="${RESUME_SHARDS:-1}"
export ARCHIVE_RUNS="${ARCHIVE_RUNS:-1}"

case "${DESIGN}" in
  gtr)
    export MODEL_PAIRS="GTR_PF06346:GTR_PF06346,GTR_PF06346_G4:GTR_PF06346_G4,GTR_PF06346_I_G4:GTR_PF06346_I_G4"
    export SCENARIO_SET=fig2
    export RUN_NAME="${RUN_NAME:-paper-matched-gtr-rate-invar-reps${REPS}-${SOURCE_ID}-${RUN_DATE}}"
    ;;
  protein)
    export MODEL_PAIRS="LG_G4:LG_G4,WAG_G4:WAG_G4,JTT_G4:JTT_G4,Q.PFAM_G4:Q.PFAM_G4"
    export SCENARIO_SET=fig2
    export RUN_NAME="${RUN_NAME:-paper-matched-protein-g4-reps${REPS}-${SOURCE_ID}-${RUN_DATE}}"
    ;;
  misspec)
    export MODEL_PAIRS="GTR_PF06346:GTR_PF06346,GTR_PF06346:JC,GTR_PF06346:K2P,GTR_PF06346:F81"
    export SCENARIO_SET=misspecification
    export RUN_NAME="${RUN_NAME:-paper-gtr-misspec-reps${REPS}-${SOURCE_ID}-${RUN_DATE}}"
    ;;
  *)
    echo "Unknown design: ${DESIGN}" >&2
    usage >&2
    exit 2
    ;;
esac

for required in RUN_NAME MODEL_PAIRS SCENARIO_SET TREE_CASES SITE_LENGTHS BRANCH_LENGTHS; do
  if [[ -z "${!required}" ]]; then
    echo "Required run-spec field is empty: ${required}" >&2
    exit 2
  fi
done

printf 'design=%s\nrun_name=%s\nreps=%s\nshards=%s\nmodels=%s\nscenarios=%s\ntrees=%s\nsites=%s\nbranches=%s\nsimulator=%s\n' \
  "${DESIGN}" "${RUN_NAME}" "${REPS}" "${SHARD_COUNT}" "${MODEL_PAIRS}" \
  "${SCENARIO_SET}" "${TREE_CASES}" "${SITE_LENGTHS}" "${BRANCH_LENGTHS}" "${SIMULATOR}"

if [[ "${DRY_RUN}" == "1" ]]; then
  exit 0
fi

sbatch --parsable "$@" "${SUBMIT_SCRIPT}"
