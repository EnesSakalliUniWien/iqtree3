#!/usr/bin/env bash
set -euo pipefail

WORK_DIR="${WORK_DIR:-/lisc/data/scratch/menche/sakalli/iq-tree-satute/full-head-to-head}"
RUN_NAME="${RUN_NAME:?Set RUN_NAME to the run root to archive}"
OUT_ROOT="${OUT_ROOT:-${WORK_DIR}/${RUN_NAME}}"
ARCHIVE_DIR="${OUT_ROOT}/run-archives"
TMP_BASE="${TMPDIR:-${ARCHIVE_DIR}}"
ZSTD_LEVEL="${ZSTD_LEVEL:-3}"
ZSTD_THREADS="${ZSTD_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"

if [[ -n "${SLURM_JOB_ID:-}" && -z "${TMPDIR:-}" ]]; then
  echo "Slurm job has no TMPDIR; refusing to stage partial archives on shared scratch" >&2
  exit 2
fi

if ! command -v zstd >/dev/null 2>&1; then
  echo "Partial-run archiving requires zstd on PATH" >&2
  exit 2
fi

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
  archive="${ARCHIVE_DIR}/runs_shard_${shard_index}_partial_job_${SLURM_JOB_ID:-manual}.tar.zst"
  checksum="${archive}.sha256"
  tmp_archive="${TMP_BASE}/.${RUN_NAME}_${shard_name}_${SLURM_JOB_ID:-manual}.tar.zst.tmp"

  if [[ -s "${archive}" ]] && zstd -q -t "${archive}" && tar -I zstd -tf "${archive}" >/dev/null 2>&1; then
    rm -rf "${runs_dir}"
    skipped=$((skipped + 1))
    removed=$((removed + 1))
    echo "Existing valid archive; removed runs: ${shard_name}"
    continue
  fi

  rm -f "${tmp_archive}" "${archive}.tmp"
  echo "Archiving ${shard_name}"
  if ! tar -C "${shard_dir}" -I "zstd -T${ZSTD_THREADS} -${ZSTD_LEVEL} -q" -cf "${tmp_archive}" runs; then
    echo "Failed to create temporary archive for ${shard_name}" >&2
    rm -f "${tmp_archive}"
    failed=$((failed + 1))
    continue
  fi
  if ! zstd -q -t "${tmp_archive}" || ! tar -I zstd -tf "${tmp_archive}" >/dev/null; then
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
  if ! zstd -q -t "${archive}" || ! tar -I zstd -tf "${archive}" >/dev/null; then
    echo "Final archive failed validation for ${shard_name}" >&2
    rm -f "${archive}" "${checksum}"
    failed=$((failed + 1))
    continue
  fi
  (
    cd "${ARCHIVE_DIR}"
    sha256sum "$(basename "${archive}")" > "$(basename "${checksum}").tmp"
    mv "$(basename "${checksum}").tmp" "$(basename "${checksum}")"
    sha256sum -c "$(basename "${checksum}")"
  )
  rm -rf "${runs_dir}"
  archived=$((archived + 1))
  removed=$((removed + 1))
done < <(find "${OUT_ROOT}/shards" -maxdepth 2 -type d -name runs -print0 2>/dev/null | sort -z)

echo "Archived new shards: ${archived}"
echo "Existing archives reused: ${skipped}"
echo "Runs directories removed: ${removed}"
echo "Failures: ${failed}"

remaining="$(find "${OUT_ROOT}/shards" -maxdepth 2 -type d -name runs 2>/dev/null | wc -l | tr -d ' ')"
archives="$(find "${ARCHIVE_DIR}" -type f \( -name 'runs_shard_*.tar.gz' -o -name 'runs_shard_*.tar.zst' \) 2>/dev/null | wc -l | tr -d ' ')"
echo "Remaining runs directories: ${remaining}"
echo "Archives present: ${archives}"

if [[ "${failed}" -ne 0 ]]; then
  exit 1
fi
