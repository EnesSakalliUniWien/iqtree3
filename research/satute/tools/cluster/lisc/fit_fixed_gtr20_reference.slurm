#!/usr/bin/env bash
#SBATCH --job-name=satute-gtr20-ref
#SBATCH --account=menche
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00

set -euo pipefail

ALIGNMENT="${ALIGNMENT:?Set ALIGNMENT to the independent protein alignment}"
TREE="${TREE:?Set TREE to the fixed reference topology}"
OUT_DIR="${OUT_DIR:?Set OUT_DIR for the fitted reference artifacts}"
IQTREE_BIN="${IQTREE_BIN:?Set IQTREE_BIN to the IQ-TREE SatuTe executable}"
LISC_RUNTIME_MODULES="${LISC_RUNTIME_MODULES:-GCCcore/14.3.0}"
MODEL="${MODEL:-GTR20+F+G4}"
SEED="${SEED:-20260713}"

if command -v module >/dev/null 2>&1 && [[ -n "${LISC_RUNTIME_MODULES}" ]]; then
  set +u
  # shellcheck disable=SC2086
  module load ${LISC_RUNTIME_MODULES}
  set -u
fi

for path in "${ALIGNMENT}" "${TREE}"; do
  if [[ ! -s "${path}" ]]; then
    echo "Missing reference input: ${path}" >&2
    exit 2
  fi
done
if [[ ! -x "${IQTREE_BIN}" ]]; then
  echo "Missing IQ-TREE executable: ${IQTREE_BIN}" >&2
  exit 2
fi

mkdir -p "${OUT_DIR}"
PREFIX="${OUT_DIR}/gtr20_reference_fit"

{
  printf "field\tvalue\n"
  printf "created_utc\t%s\n" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf "model_request\t%s\n" "${MODEL}"
  printf "seed\t%s\n" "${SEED}"
  printf "runtime_modules\t%s\n" "${LISC_RUNTIME_MODULES}"
  printf "alignment_path\t%s\n" "${ALIGNMENT}"
  printf "alignment_sha256\t%s\n" "$(sha256sum "${ALIGNMENT}" | awk '{print $1}')"
  printf "tree_path\t%s\n" "${TREE}"
  printf "tree_sha256\t%s\n" "$(sha256sum "${TREE}" | awk '{print $1}')"
  printf "iqtree_sha256\t%s\n" "$(sha256sum "${IQTREE_BIN}" | awk '{print $1}')"
} > "${OUT_DIR}/fit_manifest.tsv"

"${IQTREE_BIN}" \
  -s "${ALIGNMENT}" \
  -te "${TREE}" \
  -m "${MODEL}" \
  --prefix "${PREFIX}" \
  --seed "${SEED}" \
  -T "${SLURM_CPUS_PER_TASK}" \
  --redo

sha256sum "${PREFIX}.iqtree" "${PREFIX}.treefile" > "${OUT_DIR}/fit_output_sha256.txt"
