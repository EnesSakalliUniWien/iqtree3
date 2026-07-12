#!/usr/bin/env bash
set -euo pipefail

WORK_DIR="${WORK_DIR:-/lisc/data/scratch/menche/sakalli/iq-tree-satute/full-head-to-head}"
RUN_NAME="${RUN_NAME:?Set RUN_NAME to the run root to archive}"
OUT_ROOT="${OUT_ROOT:-${WORK_DIR}/${RUN_NAME}}"
ARCHIVE_DIR="${OUT_ROOT}/run-archives"
TMP_BASE="${TMPDIR:-${ARCHIVE_DIR}}"

if [[ ! -d "${OUT_ROOT}/shards" ]]; then
  echo "Missing shard directory: ${OUT_ROOT}/shards" >&2
  exit 2
fi

mkdir -p "${ARCHIVE_DIR}"

archived=0
removed=0
failed=0
skipped=0

while IFS= read -r -d '' runs_dir; do
  shard_dir="$(dirname "${runs_dir}")"
  shard_name="$(basename "${shard_dir}")"
  shard_index="${shard_name#shard-}"
  archive="${ARCHIVE_DIR}/runs_shard_${shard_index}_partial_job_${SLURM_JOB_ID:-manual}.tar.gz"
  tmp_archive="${TMP_BASE}/.${RUN_NAME}_${shard_name}_${SLURM_JOB_ID:-manual}.tar.gz.tmp"

  if [[ -s "${archive}" ]] && tar -tzf "${archive}" >/dev/null 2>&1; then
    rm -rf "${runs_dir}"
    skipped=$((skipped + 1))
    removed=$((removed + 1))
    echo "Existing valid archive; removed runs: ${shard_name}"
    continue
  fi

  rm -f "${tmp_archive}" "${archive}.tmp"
  echo "Archiving ${shard_name}"
  if ! tar -C "${shard_dir}" -czf "${tmp_archive}" runs; then
    echo "Failed to create temporary archive for ${shard_name}" >&2
    rm -f "${tmp_archive}"
    failed=$((failed + 1))
    continue
  fi
  if ! tar -tzf "${tmp_archive}" >/dev/null; then
    echo "Temporary archive failed validation for ${shard_name}" >&2
    rm -f "${tmp_archive}"
    failed=$((failed + 1))
    continue
  fi
  if ! cp "${tmp_archive}" "${archive}.tmp"; then
    echo "Failed to copy archive to ${ARCHIVE_DIR} for ${shard_name}" >&2
    rm -f "${tmp_archive}" "${archive}.tmp"
    failed=$((failed + 1))
    continue
  fi
  mv "${archive}.tmp" "${archive}"
  rm -f "${tmp_archive}"
  if ! tar -tzf "${archive}" >/dev/null; then
    echo "Final archive failed validation for ${shard_name}" >&2
    failed=$((failed + 1))
    continue
  fi
  rm -rf "${runs_dir}"
  archived=$((archived + 1))
  removed=$((removed + 1))
done < <(find "${OUT_ROOT}/shards" -maxdepth 2 -type d -name runs -print0 2>/dev/null | sort -z)

echo "Archived new shards: ${archived}"
echo "Existing archives reused: ${skipped}"
echo "Runs directories removed: ${removed}"
echo "Failures: ${failed}"

remaining="$(find "${OUT_ROOT}/shards" -maxdepth 2 -type d -name runs 2>/dev/null | wc -l | tr -d ' ')"
archives="$(find "${ARCHIVE_DIR}" -name 'runs_shard_*.tar.gz' -type f 2>/dev/null | wc -l | tr -d ' ')"
echo "Remaining runs directories: ${remaining}"
echo "Archives present: ${archives}"

if [[ "${failed}" -ne 0 ]]; then
  exit 1
fi
