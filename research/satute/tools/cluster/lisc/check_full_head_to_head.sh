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
find "${OUT_ROOT}/shards" -path '*/head_to_head_detail.tsv' -type f 2>/dev/null | wc -l || true

echo
echo "Shard archives:"
find "${OUT_ROOT}/run-archives" -type f \
  \( -name 'runs_shard_*.tar.gz' -o -name 'runs_shard_*.tar.zst' \) \
  2>/dev/null | wc -l || true

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
