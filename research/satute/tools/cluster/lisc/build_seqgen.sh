#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
ROOT_DIR="${ROOT_DIR:-$(cd "${PROJECT_DIR}/../.." && pwd)}"
WORK_DIR="${WORK_DIR:-/lisc/data/scratch/menche/sakalli/iq-tree-satute/full-head-to-head}"
SEQGEN_SRC="${SEQGEN_SRC:-}"
SEQGEN_BUILD_DIR="${SEQGEN_BUILD_DIR:-${WORK_DIR}/seq-gen-1.3.4-source}"

if [[ -z "${SEQGEN_SRC}" ]]; then
  for candidate in \
    "${ROOT_DIR}/tools/seq-gen-1.3.4/source" \
    "/lisc/home/user/sakalli/projects/satute-eigenvector-reanalysis/tools/seq-gen-1.3.4/source"; do
    if [[ -f "${candidate}/seq-gen.c" && -f "${candidate}/Makefile" ]]; then
      SEQGEN_SRC="${candidate}"
      break
    fi
  done
fi

if [[ -z "${SEQGEN_SRC}" || ! -f "${SEQGEN_SRC}/Makefile" ]]; then
  echo "Cannot find Seq-Gen source. Set SEQGEN_SRC to the seq-gen-1.3.4/source directory." >&2
  exit 2
fi

mkdir -p "${WORK_DIR}"
if [[ "${SEQGEN_SRC}" != "${SEQGEN_BUILD_DIR}" ]]; then
  mkdir -p "${SEQGEN_BUILD_DIR}"
  rsync -a "${SEQGEN_SRC}/" "${SEQGEN_BUILD_DIR}/"
fi

make -C "${SEQGEN_BUILD_DIR}"
printf "%s\n" "${SEQGEN_BUILD_DIR}/seq-gen"
