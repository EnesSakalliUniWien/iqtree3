#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
ROOT_DIR="${ROOT_DIR:-$(cd "${PROJECT_DIR}/../.." && pwd)}"
WORK_DIR="${WORK_DIR:-/lisc/data/scratch/menche/sakalli/iq-tree-satute/full-head-to-head}"
RUN_NAME="${RUN_NAME:?RUN_NAME is required}"
OUT_ROOT="${OUT_ROOT:-${WORK_DIR}/${RUN_NAME}}"
IQTREE_BIN="${IQTREE_BIN:-${WORK_DIR}/iqtree3-build/iqtree3}"
EVONAPS_BRANCH_LENGTHS="${EVONAPS_BRANCH_LENGTHS:-${PROJECT_DIR}/references/evonaps/evonaps_16taxon_branch_lengths.tsv}"
LISC_RUNTIME_MODULES="${LISC_RUNTIME_MODULES:-SciPy-bundle/2025.06-gfbf-2025a}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
DESIGN="${DESIGN:?DESIGN is required}"
REPS="${REPS:?REPS is required}"
TRAINING_REPS="${TRAINING_REPS:-0}"
MODELS="${MODELS:-JC,GTR_PF06346}"
SITE_LENGTHS="${SITE_LENGTHS:-100,1000}"
BRANCH_LENGTHS="${BRANCH_LENGTHS:-4,8}"
TREE_CASES="${TREE_CASES:-five_external,sixteen_internal}"
SHARD_INDEX="${SHARD_INDEX:-${SLURM_ARRAY_TASK_ID:-0}}"
SHARD_COUNT="${SHARD_COUNT:-1000}"
ZSTD_LEVEL="${ZSTD_LEVEL:-3}"
ZSTD_THREADS="${ZSTD_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"

if command -v module >/dev/null 2>&1; then
  set +u
  # shellcheck disable=SC2086
  module load ${LISC_RUNTIME_MODULES}
  set -u
fi
if [[ ! -x "${IQTREE_BIN}" ]]; then
  echo "Missing IQ-TREE binary: ${IQTREE_BIN}" >&2
  exit 2
fi
if [[ ! -s "${EVONAPS_BRANCH_LENGTHS}" ]]; then
  echo "Missing EvoNAPS table: ${EVONAPS_BRANCH_LENGTHS}" >&2
  exit 2
fi
if [[ -n "${SLURM_JOB_ID:-}" && -z "${TMPDIR:-}" ]]; then
  echo "Slurm job has no TMPDIR" >&2
  exit 2
fi

SHARD_DIR="${OUT_ROOT}/shards/shard-${SHARD_INDEX}"
WORK_SHARD="${TMPDIR:-${SHARD_DIR}}/satute-fdr-${SLURM_ARRAY_JOB_ID:-manual}-${SHARD_INDEX}"
mkdir -p "${OUT_ROOT}/run-archives" "${SHARD_DIR}" "${WORK_SHARD}"

copy_atomic() {
  local source="$1" destination="$2" temporary="${destination}.tmp.${SLURM_JOB_ID:-$$}"
  cp "${source}" "${temporary}"
  mv "${temporary}" "${destination}"
}

cleanup() {
  local status=$?
  trap - EXIT
  if [[ "${status}" -ne 0 ]]; then
    for name in truth_branch_audit.tsv replicate_fdp.tsv calibration_summary.tsv timing.tsv; do
      if [[ -s "${WORK_SHARD}/${name}" ]]; then
        copy_atomic "${WORK_SHARD}/${name}" "${SHARD_DIR}/${name}"
      fi
    done
  fi
  if [[ -n "${TMPDIR:-}" ]]; then
    rm -rf "${WORK_SHARD}"
  fi
  exit "${status}"
}
trap cleanup EXIT

if [[ "${SHARD_INDEX}" == "0" ]]; then
  matched_model_pairs="$(
    tr ',' '\n' <<<"${MODELS}" | awk 'NF {gsub(/^[[:space:]]+|[[:space:]]+$/, ""); printf "%s%s:%s", separator, $0, $0; separator=","}'
  )"
  resolved_model_pairs="$(
    PYTHONPATH="${PROJECT_DIR}/src${PYTHONPATH:+:${PYTHONPATH}}" "${PYTHON_BIN}" \
      "${PROJECT_DIR}/experiments/002_relative_weighting/run.py" \
      --model-pairs "${matched_model_pairs}" \
      --print-resolved-model-pairs
  )"
  manifest_tmp="${OUT_ROOT}/.manifest.${SLURM_JOB_ID:-manual}.tmp"
  {
    printf 'field\tvalue\n'
    printf 'created_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'schema_version\t2\n'
    printf 'design\t%s\n' "${DESIGN}"
    printf 'replicates\t%s\n' "${REPS}"
    printf 'training_replicates\t%s\n' "${TRAINING_REPS}"
    printf 'shard_count\t%s\n' "${SHARD_COUNT}"
    printf 'models\t%s\n' "${MODELS}"
    printf 'resolved_model_pairs\t%s\n' "${resolved_model_pairs}"
    printf 'site_lengths\t%s\n' "${SITE_LENGTHS}"
    printf 'branch_lengths\t%s\n' "${BRANCH_LENGTHS}"
    printf 'tree_cases\t%s\n' "${TREE_CASES}"
    printf 'truth_key\tnormalized_split\n'
    printf 'runtime_modules\t%s\n' "${LISC_RUNTIME_MODULES}"
    printf 'iqtree_sha256\t%s\n' "$(sha256sum "${IQTREE_BIN}" | awk '{print $1}')"
    printf 'driver_sha256\t%s\n' "$(sha256sum "${PROJECT_DIR}/experiments/002_relative_weighting/fdr_calibration_run.py" | awk '{print $1}')"
    printf 'analysis_module_sha256\t%s\n' "$(sha256sum "${PROJECT_DIR}/src/satute_analysis/fdr_calibration.py" | awk '{print $1}')"
    printf 'config_sha256\t%s\n' "$(sha256sum "${PROJECT_DIR}/config/experiments/002_relative_weighting.yaml" | awk '{print $1}')"
  } > "${manifest_tmp}"
  mv "${manifest_tmp}" "${OUT_ROOT}/run_manifest.tsv"
fi

PYTHONPATH="${PROJECT_DIR}/src${PYTHONPATH:+:${PYTHONPATH}}" "${PYTHON_BIN}" \
  "${PROJECT_DIR}/experiments/002_relative_weighting/fdr_calibration_run.py" \
  --iqtree "${IQTREE_BIN}" \
  --outdir "${WORK_SHARD}" \
  --design "${DESIGN}" \
  --reps "${REPS}" \
  --training-reps "${TRAINING_REPS}" \
  --models "${MODELS}" \
  --site-lengths "${SITE_LENGTHS}" \
  --branch-lengths "${BRANCH_LENGTHS}" \
  --tree-cases "${TREE_CASES}" \
  --evonaps-branch-lengths "${EVONAPS_BRANCH_LENGTHS}" \
  --shard-index "${SHARD_INDEX}" \
  --shard-count "${SHARD_COUNT}"

archive_name="runs_shard_${SHARD_INDEX}_job_${SLURM_JOB_ID:-manual}.tar.zst"
archive_local="${WORK_SHARD}/${archive_name}"
tar -C "${WORK_SHARD}" -I "zstd -T${ZSTD_THREADS} -${ZSTD_LEVEL} -q" -cf "${archive_local}" runs
zstd -q -t "${archive_local}"
tar -I zstd -tf "${archive_local}" >/dev/null
(
  cd "${WORK_SHARD}"
  sha256sum "${archive_name}" > "${archive_name}.sha256"
)
copy_atomic "${archive_local}" "${OUT_ROOT}/run-archives/${archive_name}"
copy_atomic "${archive_local}.sha256" "${OUT_ROOT}/run-archives/${archive_name}.sha256"

for name in truth_branch_audit.tsv replicate_fdp.tsv calibration_summary.tsv timing.tsv; do
  test -s "${WORK_SHARD}/${name}"
  copy_atomic "${WORK_SHARD}/${name}" "${SHARD_DIR}/${name}"
done
rm -rf "${WORK_SHARD}/runs"
