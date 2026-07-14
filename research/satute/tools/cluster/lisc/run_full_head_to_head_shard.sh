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
LISC_RUNTIME_MODULES="${LISC_RUNTIME_MODULES:-GCCcore/14.3.0}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
REPS="${REPS:-1000}"
SIMULATION_MODELS="${SIMULATION_MODELS:-JC}"
EVALUATION_MODELS="${EVALUATION_MODELS:-${SIMULATION_MODELS}}"
MODEL_PAIRS="${MODEL_PAIRS:-}"
SCENARIO_SET="${SCENARIO_SET:-fig2}"
TREE_CASES="${TREE_CASES:-five_external,sixteen_internal}"
SHARD_INDEX="${SHARD_INDEX:-${SLURM_ARRAY_TASK_ID:-0}}"
SHARD_COUNT="${SHARD_COUNT:-1000}"
ARCHIVE_RUNS="${ARCHIVE_RUNS:-1}"
RESUME_SHARDS="${RESUME_SHARDS:-1}"
ZSTD_LEVEL="${ZSTD_LEVEL:-3}"
ZSTD_THREADS="${ZSTD_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
SIMULATOR="${SIMULATOR:-seq-gen}"
PAPER_GRID="${PAPER_GRID:-1}"
SITE_LENGTHS="${SITE_LENGTHS:-100}"
BRANCH_LENGTHS="${BRANCH_LENGTHS:-0.1,1.0}"

for boolean_name in ARCHIVE_RUNS RESUME_SHARDS; do
  boolean_value="${!boolean_name}"
  case "${boolean_value}" in
    0|1|false|true) ;;
    *)
      echo "${boolean_name} must be 0, 1, false, or true; got: ${boolean_value}" >&2
      exit 2
      ;;
  esac
done

if command -v module >/dev/null 2>&1 && [[ -n "${LISC_RUNTIME_MODULES}" ]]; then
  set +u
  # shellcheck disable=SC2086
  module load ${LISC_RUNTIME_MODULES}
  set -u
fi

if [[ ! -x "${IQTREE_BIN}" ]]; then
  echo "Missing IQ-TREE binary: ${IQTREE_BIN}" >&2
  exit 2
fi

if [[ "${SIMULATOR}" == "seq-gen" && ! -x "${SEQGEN_BIN}" ]]; then
  echo "Missing Seq-Gen binary: ${SEQGEN_BIN}" >&2
  exit 2
fi

if [[ "${ARCHIVE_RUNS}" == "1" || "${ARCHIVE_RUNS}" == "true" ]]; then
  if ! command -v zstd >/dev/null 2>&1; then
    echo "ARCHIVE_RUNS requires zstd on PATH" >&2
    exit 2
  fi
fi

if [[ ! -s "${EVONAPS_BRANCH_LENGTHS}" ]]; then
  echo "Missing EvoNAPS branch-length table: ${EVONAPS_BRANCH_LENGTHS}" >&2
  exit 2
fi

SHARD_DIR="${OUT_ROOT}/shards/shard-${SHARD_INDEX}"
mkdir -p "${OUT_ROOT}/logs" "${OUT_ROOT}/run-archives" "${SHARD_DIR}"

if [[ -n "${SLURM_JOB_ID:-}" && -z "${TMPDIR:-}" ]]; then
  echo "Slurm job has no TMPDIR; refusing to create IQ-TREE intermediates on shared scratch" >&2
  exit 2
fi

if [[ -n "${TMPDIR:-}" ]]; then
  WORK_SHARD_DIR="${TMPDIR}/satute-${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}-${SHARD_INDEX}"
  NODE_LOCAL_WORK=1
else
  WORK_SHARD_DIR="${SHARD_DIR}"
  NODE_LOCAL_WORK=0
fi
mkdir -p "${WORK_SHARD_DIR}"

copy_atomic() {
  local source="$1"
  local destination="$2"
  local temporary="${destination}.tmp.${SLURM_JOB_ID:-$$}"
  rm -f "${temporary}"
  cp "${source}" "${temporary}"
  mv "${temporary}" "${destination}"
}

persist_tabular_results() {
  local name
  for name in head_to_head_detail.tsv head_to_head_summary.tsv; do
    if [[ -s "${WORK_SHARD_DIR}/${name}" ]]; then
      copy_atomic "${WORK_SHARD_DIR}/${name}" "${SHARD_DIR}/${name}"
    fi
  done
}

cleanup_workdir() {
  local status=$?
  trap - EXIT
  if [[ "${status}" -ne 0 ]]; then
    set +e
    if [[ -n "${DURABLE_ARCHIVE_IN_PROGRESS:-}" ]]; then
      rm -f "${DURABLE_ARCHIVE_IN_PROGRESS}" "${DURABLE_ARCHIVE_IN_PROGRESS}.sha256"
    fi
    if [[ ! -s "${WORK_SHARD_DIR}/head_to_head_summary.tsv" || "${ARCHIVE_RUNS}" == "0" || "${ARCHIVE_RUNS}" == "false" ]]; then
      persist_tabular_results
    else
      echo "Completed table was not persisted because its required runs archive failed" >&2
    fi
    set -e
  fi
  if [[ "${NODE_LOCAL_WORK}" == "1" ]]; then
    rm -rf "${WORK_SHARD_DIR}"
  fi
  exit "${status}"
}
trap cleanup_workdir EXIT

COPY_RESUME_STATE="${RESUME_SHARDS}"
if [[ ( "${COPY_RESUME_STATE}" == "1" || "${COPY_RESUME_STATE}" == "true" ) \
      && ( "${ARCHIVE_RUNS}" == "1" || "${ARCHIVE_RUNS}" == "true" ) \
      && -s "${SHARD_DIR}/head_to_head_summary.tsv" ]] \
   && ! compgen -G "${OUT_ROOT}/run-archives/runs_shard_${SHARD_INDEX}_job_*.tar.*" >/dev/null; then
  COPY_RESUME_STATE=0
  echo "Complete table has no runs archive; forcing a full shard rerun"
fi

if [[ "${NODE_LOCAL_WORK}" == "1" \
      && ( "${COPY_RESUME_STATE}" == "1" || "${COPY_RESUME_STATE}" == "true" ) ]]; then
  for name in head_to_head_detail.tsv head_to_head_summary.tsv; do
    if [[ -s "${SHARD_DIR}/${name}" ]]; then
      cp "${SHARD_DIR}/${name}" "${WORK_SHARD_DIR}/${name}"
    fi
  done
fi

if [[ "${SHARD_INDEX}" == "0" ]]; then
  MANIFEST_TMP="${OUT_ROOT}/.run_manifest_${SLURM_JOB_ID:-manual}.tmp"
  {
    printf "field\tvalue\n"
    printf "created_utc\t%s\n" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf "slurm_array_job_id\t%s\n" "${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
    printf "run_name\t%s\n" "${RUN_NAME}"
    printf "replicates\t%s\n" "${REPS}"
    printf "shard_count\t%s\n" "${SHARD_COUNT}"
    printf "resume_shards\t%s\n" "${RESUME_SHARDS}"
    printf "tree_cases\t%s\n" "${TREE_CASES}"
    printf "simulation_models\t%s\n" "${SIMULATION_MODELS}"
    printf "evaluation_models\t%s\n" "${EVALUATION_MODELS}"
    printf "model_pairs\t%s\n" "${MODEL_PAIRS}"
    printf "scenario_set\t%s\n" "${SCENARIO_SET}"
    printf "simulator\t%s\n" "${SIMULATOR}"
    printf "paper_grid\t%s\n" "${PAPER_GRID}"
    printf "site_lengths\t%s\n" "${SITE_LENGTHS}"
    printf "branch_lengths\t%s\n" "${BRANCH_LENGTHS}"
    printf "runtime_modules\t%s\n" "${LISC_RUNTIME_MODULES}"
    printf "iqtree_sha256\t%s\n" "$(sha256sum "${IQTREE_BIN}" | awk '{print $1}')"
    if [[ "${SIMULATOR}" == "seq-gen" ]]; then
      printf "seqgen_sha256\t%s\n" "$(sha256sum "${SEQGEN_BIN}" | awk '{print $1}')"
    else
      printf "seqgen_sha256\tnot_used\n"
    fi
    printf "evonaps_sha256\t%s\n" "$(sha256sum "${EVONAPS_BRANCH_LENGTHS}" | awk '{print $1}')"
    printf "driver_sha256\t%s\n" "$(sha256sum "${PROJECT_DIR}/experiments/002_relative_weighting/run.py" | awk '{print $1}')"
    printf "satute_facade_sha256\t%s\n" "$(sha256sum "${ROOT_DIR}/tree/satute.cpp" | awk '{print $1}')"
    printf "satute_statistics_sha256\t%s\n" "$(sha256sum "${ROOT_DIR}/tree/satute/satute_statistics.cpp" | awk '{print $1}')"
    printf "satute_support_sha256\t%s\n" "$(sha256sum "${ROOT_DIR}/tree/satute/satute_support.cpp" | awk '{print $1}')"
    printf "satute_report_writer_sha256\t%s\n" "$(sha256sum "${ROOT_DIR}/tree/satute/satute_report_writer.cpp" | awk '{print $1}')"
  } > "${MANIFEST_TMP}"
  mv "${MANIFEST_TMP}" "${OUT_ROOT}/run_manifest.tsv"
fi

model_pair_args=()
if [[ -n "${MODEL_PAIRS}" ]]; then
  model_pair_args=(--model-pairs "${MODEL_PAIRS}")
fi

grid_args=()
if [[ "${PAPER_GRID}" == "1" || "${PAPER_GRID}" == "true" ]]; then
  grid_args=(--paper-grid)
else
  grid_args=(--site-lengths "${SITE_LENGTHS}" --branch-lengths "${BRANCH_LENGTHS}")
fi

PYTHONPATH="${PROJECT_DIR}/src${PYTHONPATH:+:${PYTHONPATH}}" "${PYTHON_BIN}" \
  "${PROJECT_DIR}/experiments/002_relative_weighting/run.py" \
  --iqtree "${IQTREE_BIN}" \
  --seqgen "${SEQGEN_BIN}" \
  --simulator "${SIMULATOR}" \
  --evonaps-branch-lengths "${EVONAPS_BRANCH_LENGTHS}" \
  --outdir "${WORK_SHARD_DIR}" \
  --reps "${REPS}" \
  "${grid_args[@]}" \
  --tree-cases "${TREE_CASES}" \
  --simulation-models "${SIMULATION_MODELS}" \
  --evaluation-models "${EVALUATION_MODELS}" \
  "${model_pair_args[@]}" \
  --scenario-set "${SCENARIO_SET}" \
  --shard-index "${SHARD_INDEX}" \
  --shard-count "${SHARD_COUNT}" \
  --resume

if [[ "${ARCHIVE_RUNS}" == "1" || "${ARCHIVE_RUNS}" == "true" ]]; then
  RUNS_DIR="${WORK_SHARD_DIR}/runs"
  if [[ -d "${RUNS_DIR}" ]]; then
    ARCHIVE_BASENAME="runs_shard_${SHARD_INDEX}_job_${SLURM_JOB_ID:-manual}.tar.zst"
    LOCAL_ARCHIVE="${WORK_SHARD_DIR}/${ARCHIVE_BASENAME}"
    LOCAL_CHECKSUM="${LOCAL_ARCHIVE}.sha256"
    ARCHIVE="${OUT_ROOT}/run-archives/${ARCHIVE_BASENAME}"
    CHECKSUM="${ARCHIVE}.sha256"
    DURABLE_ARCHIVE_IN_PROGRESS="${ARCHIVE}"

    tar -C "${WORK_SHARD_DIR}" -I "zstd -T${ZSTD_THREADS} -${ZSTD_LEVEL} -q" -cf "${LOCAL_ARCHIVE}" runs
    zstd -q -t "${LOCAL_ARCHIVE}"
    tar -I zstd -tf "${LOCAL_ARCHIVE}" >/dev/null
    (
      cd "${WORK_SHARD_DIR}"
      sha256sum "${ARCHIVE_BASENAME}" > "${ARCHIVE_BASENAME}.sha256"
    )

    copy_atomic "${LOCAL_ARCHIVE}" "${ARCHIVE}"
    copy_atomic "${LOCAL_CHECKSUM}" "${CHECKSUM}"
    zstd -q -t "${ARCHIVE}"
    tar -I zstd -tf "${ARCHIVE}" >/dev/null
    (
      cd "${OUT_ROOT}/run-archives"
      sha256sum -c "${ARCHIVE_BASENAME}.sha256"
    )
    DURABLE_ARCHIVE_IN_PROGRESS=""
    rm -rf "${RUNS_DIR}"
    echo "Archived and verified shard runs: ${ARCHIVE}"
  fi
fi

persist_tabular_results
