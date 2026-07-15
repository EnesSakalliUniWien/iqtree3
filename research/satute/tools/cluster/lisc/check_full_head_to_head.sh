#!/usr/bin/env bash
set -euo pipefail

WORK_DIR="${WORK_DIR:-/lisc/data/scratch/menche/sakalli/iq-tree-satute/full-head-to-head}"
RUN_NAME="${RUN_NAME:-paper-fig2-jc-reps1000}"
OUT_ROOT="${OUT_ROOT:-${WORK_DIR}/${RUN_NAME}}"
SHARD_COUNT="${SHARD_COUNT:-1000}"

echo "Queue:"
squeue --me || true

echo
echo "Shard detail files:"
detail_count=$(find "${OUT_ROOT}/shards" -path '*/head_to_head_detail.tsv' -type f 2>/dev/null | wc -l | tr -d ' ')
echo "${detail_count}"

for companion in head_to_head_summary.tsv head_to_head_branch_audit.tsv head_to_head_fdr_audit.tsv head_to_head_timing.tsv; do
  count=$(find "${OUT_ROOT}/shards" -path "*/${companion}" -type f 2>/dev/null | wc -l | tr -d ' ')
  echo "${companion}: ${count}"
  if [[ "${count}" -ne "${SHARD_COUNT}" ]]; then
    echo "ERROR: expected ${SHARD_COUNT} ${companion} files" >&2
    exit 1
  fi
done
if [[ "${detail_count}" -ne "${SHARD_COUNT}" ]]; then
  echo "ERROR: expected ${SHARD_COUNT} detail files" >&2
  exit 1
fi

echo
echo "Shard archives:"
archive_count=$(find "${OUT_ROOT}/run-archives" -type f \
  \( -name 'runs_shard_*.tar.gz' -o -name 'runs_shard_*.tar.zst' \) \
  2>/dev/null | wc -l | tr -d ' ')
echo "${archive_count}"
checksum_count=$(find "${OUT_ROOT}/run-archives" -type f -name 'runs_shard_*.tar.*.sha256' 2>/dev/null | wc -l | tr -d ' ')
echo "Archive checksums: ${checksum_count}"
if [[ "${archive_count}" -ne "${SHARD_COUNT}" || "${checksum_count}" -ne "${SHARD_COUNT}" ]]; then
  echo "ERROR: archive/checksum counts do not match expected shards" >&2
  exit 1
fi
runs_dirs=$(find "${OUT_ROOT}" -type d -name runs 2>/dev/null | wc -l | tr -d ' ')
echo "Durable runs directories: ${runs_dirs}"
if [[ "${runs_dirs}" -ne 0 ]]; then
  echo "ERROR: durable runs directories remain under ${OUT_ROOT}" >&2
  exit 1
fi

if [[ "${archive_count}" -gt 0 ]]; then
  sample_archive=$(find "${OUT_ROOT}/run-archives" -type f -name 'runs_shard_*.tar.zst' | sort | head -1)
  if [[ -n "${sample_archive}" ]]; then
    command -v zstd >/dev/null 2>&1 || { echo "ERROR: zstd is required to validate archives" >&2; exit 1; }
    zstd -t "${sample_archive}"
    tar -I zstd -tf "${sample_archive}" >/dev/null
    checksum_file="${sample_archive}.sha256"
    (cd "$(dirname "${sample_archive}")" && sha256sum -c "$(basename "${checksum_file}")")
  fi
fi

echo
echo "Expected shards:"
echo "${SHARD_COUNT}"

echo
echo "Detail lines including headers:"
find "${OUT_ROOT}/shards" -path '*/head_to_head_detail.tsv' -type f -print0 2>/dev/null \
  | xargs -0 wc -l 2>/dev/null || true

echo
echo "Scratch usage:"
du -sh "${OUT_ROOT}" 2>/dev/null || true
