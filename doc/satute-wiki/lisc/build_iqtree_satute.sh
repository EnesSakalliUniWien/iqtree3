#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
WORK_DIR="${WORK_DIR:-/lisc/data/scratch/menche/sakalli/iq-tree-satute/full-head-to-head}"
BUILD_DIR="${BUILD_DIR:-${WORK_DIR}/iqtree3-build}"
BUILD_THREADS="${BUILD_THREADS:-8}"
LISC_BUILD_MODULES="${LISC_BUILD_MODULES:-CMake/3.31.8-GCCcore-14.3.0 Eigen/3.4.0-GCCcore-14.3.0 Boost/1.88.0-GCC-14.3.0}"

if command -v module >/dev/null 2>&1 && [[ -n "${LISC_BUILD_MODULES}" ]]; then
  set +u
  # shellcheck disable=SC2086
  module load ${LISC_BUILD_MODULES}
  set -u
fi

mkdir -p "${WORK_DIR}" "${BUILD_DIR}"
cmake -S "${ROOT_DIR}" -B "${BUILD_DIR}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DUSE_CMAPLE=OFF \
  -DUSE_CMAPLE_AA=OFF \
  -DUSE_LSD2=OFF
cmake --build "${BUILD_DIR}" --target iqtree3 -j "${BUILD_THREADS}"

printf "%s\n" "${BUILD_DIR}/iqtree3"
