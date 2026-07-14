#!/usr/bin/env bash
#SBATCH --job-name=satute-clean-failed
#SBATCH --account=menche
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=04:00:00
#SBATCH --output=research/satute/artifacts/cluster/lisc/logs/cleanup_%j.out
#SBATCH --error=research/satute/artifacts/cluster/lisc/logs/cleanup_%j.err

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(cd "${SLURM_SUBMIT_DIR:-.}/research/satute" && pwd)}"
WORK_DIR="${WORK_DIR:-/lisc/data/scratch/menche/sakalli/iq-tree-satute/full-head-to-head}"
MODE="${MODE:-dry-run}"
PARALLELISM="${PARALLELISM:-${SLURM_CPUS_PER_TASK:-4}}"
REPORT_DIR="${REPORT_DIR:-${PROJECT_DIR}/artifacts/cluster/lisc/recovery-${SLURM_JOB_ID:-manual}}"
MANIFEST="${REPORT_DIR}/failed_shard_runs.tsv"
MANIFEST_TMP="${MANIFEST}.tmp"

case "${MODE}" in
  dry-run|delete) ;;
  *)
    echo "MODE must be dry-run or delete, got: ${MODE}" >&2
    exit 2
    ;;
esac

if ! [[ "${PARALLELISM}" =~ ^[1-9][0-9]*$ ]]; then
  echo "PARALLELISM must be a positive integer, got: ${PARALLELISM}" >&2
  exit 2
fi

mkdir -p "${PROJECT_DIR}/artifacts/cluster/lisc/logs" "${REPORT_DIR}"
printf "array_job_id\trun_name\tshard_index\truns_dir\tdetail_path\tdetail_bytes\n" > "${MANIFEST_TMP}"

collect_failed_shards() {
  local array_job_id="$1"
  local run_name="$2"
  local expected_failed="$3"
  local failed_count=0
  local runs_count=0
  local detail_count=0

  while IFS='|' read -r job_id state; do
    case "${job_id}" in
      "${array_job_id}"_*) ;;
      *) continue ;;
    esac

    local shard_index="${job_id#${array_job_id}_}"
    local shard_dir="${WORK_DIR}/${run_name}/shards/shard-${shard_index}"
    local runs_dir="${shard_dir}/runs"
    local detail_path="${shard_dir}/head_to_head_detail.tsv"
    local detail_bytes=0

    failed_count=$((failed_count + 1))
    if [[ -d "${runs_dir}" ]]; then
      runs_count=$((runs_count + 1))
    fi
    if [[ -s "${detail_path}" ]]; then
      detail_count=$((detail_count + 1))
      detail_bytes="$(stat -c '%s' "${detail_path}")"
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
      "${array_job_id}" "${run_name}" "${shard_index}" \
      "${runs_dir}" "${detail_path}" "${detail_bytes}" >> "${MANIFEST_TMP}"
  done < <(sacct -X --array -j "${array_job_id}" -s FAILED -n -P -o JobID,State)

  printf "Validated job %s: failed=%s runs=%s nonempty_detail=%s expected=%s\n" \
    "${array_job_id}" "${failed_count}" "${runs_count}" "${detail_count}" "${expected_failed}"

  if [[ "${failed_count}" -ne "${expected_failed}" ]]; then
    echo "Unexpected failed-task count for ${array_job_id}" >&2
    exit 1
  fi
  if [[ "${runs_count}" -ne "${expected_failed}" ]]; then
    echo "Not every failed task has exactly one runs directory for ${array_job_id}" >&2
    exit 1
  fi
  if [[ "${detail_count}" -ne "${expected_failed}" ]]; then
    echo "Not every failed task has a recoverable detail file for ${array_job_id}" >&2
    exit 1
  fi
}

collect_failed_shards 5072503 full-lg-native-reps1000-20260713 801
collect_failed_shards 5072504 full-gtr-pf06346-misspec-native-reps1000-20260713 40

mv "${MANIFEST_TMP}" "${MANIFEST}"
sha256sum "${MANIFEST}" > "${MANIFEST}.sha256"

manifest_rows="$(( $(wc -l < "${MANIFEST}") - 1 ))"
if [[ "${manifest_rows}" -ne 841 ]]; then
  echo "Cleanup manifest must contain exactly 841 rows, got ${manifest_rows}" >&2
  exit 1
fi

echo "Cleanup manifest: ${MANIFEST}"
echo "Cleanup manifest rows: ${manifest_rows}"
echo "Mode: ${MODE}"

if [[ "${MODE}" == "dry-run" ]]; then
  echo "Dry run complete; no files were removed."
  exit 0
fi

awk -F '\t' 'NR > 1 { printf "%s%c", $4, 0 }' "${MANIFEST}" \
  | xargs -0 -r -n 1 -P "${PARALLELISM}" rm -rf --

remaining=0
while IFS=$'\t' read -r array_job_id run_name shard_index runs_dir detail_path detail_bytes; do
  [[ "${array_job_id}" == "array_job_id" ]] && continue
  if [[ -e "${runs_dir}" ]]; then
    echo "Runs directory still exists: ${runs_dir}" >&2
    remaining=$((remaining + 1))
  fi
  if [[ ! -s "${detail_path}" ]]; then
    echo "Recoverable detail file was lost or emptied: ${detail_path}" >&2
    exit 1
  fi
done < "${MANIFEST}"

if [[ "${remaining}" -ne 0 ]]; then
  echo "Cleanup left ${remaining} runs directories" >&2
  exit 1
fi

echo "Deleted failed-shard runs directories: ${manifest_rows}"
/usr/lpp/mmfs/bin/mmlsquota -u "${USER}" data
