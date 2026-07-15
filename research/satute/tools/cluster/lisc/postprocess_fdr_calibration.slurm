#!/usr/bin/env bash
#SBATCH --job-name=satute-fdr-post
#SBATCH --account=menche
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=04:00:00
#SBATCH --output=research/satute/artifacts/cluster/lisc/logs/fdr_post_%j.out
#SBATCH --error=research/satute/artifacts/cluster/lisc/logs/fdr_post_%j.err

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "${SLURM_SUBMIT_DIR:-.}" && pwd)}"
PROJECT_DIR="${PROJECT_DIR:-${ROOT_DIR}/research/satute}"
WORK_DIR="${WORK_DIR:-/lisc/data/scratch/menche/sakalli/iq-tree-satute/full-head-to-head}"
RUN_ROOT="${RUN_ROOT:?RUN_ROOT is required}"
OUTPUT_ROOT="${OUTPUT_ROOT:?OUTPUT_ROOT is required}"
ARRAY_JOB_ID="${ARRAY_JOB_ID:?ARRAY_JOB_ID is required}"
EXPECTED_REPS="${EXPECTED_REPS:?EXPECTED_REPS is required}"
TRAINING_REPS="${TRAINING_REPS:-0}"
EXPECTED_SHARDS="${EXPECTED_SHARDS:-1000}"
SOURCE_ID="${SOURCE_ID:-unknown}"
LISC_RUNTIME_MODULES="${LISC_RUNTIME_MODULES:-SciPy-bundle/2025.06-gfbf-2025a}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
ZSTD_LEVEL="${ZSTD_LEVEL:-7}"
ZSTD_THREADS="${ZSTD_THREADS:-${SLURM_CPUS_PER_TASK:-2}}"
LOG_DIR="${PROJECT_DIR}/artifacts/cluster/lisc/logs"
INPUT_LOG_DIR="${INPUT_LOG_DIR:-${LOG_DIR}}"

if command -v module >/dev/null 2>&1; then
  set +u
  # shellcheck disable=SC2086
  module load ${LISC_RUNTIME_MODULES}
  set -u
fi
mkdir -p "${LOG_DIR}"
if [[ ! -d "${RUN_ROOT}/shards" || ! -s "${RUN_ROOT}/run_manifest.tsv" ]]; then
  echo "Incomplete calibration run root: ${RUN_ROOT}" >&2
  exit 2
fi
if [[ -e "${OUTPUT_ROOT}" ]]; then
  echo "Refusing existing postprocess output: ${OUTPUT_ROOT}" >&2
  exit 2
fi

count_nonempty() {
  find "$1" -type f -name "$2" -size +0c | wc -l
}

for table in truth_branch_audit.tsv replicate_fdp.tsv calibration_summary.tsv timing.tsv; do
  count="$(count_nonempty "${RUN_ROOT}/shards" "${table}")"
  if [[ "${count}" -ne "${EXPECTED_SHARDS}" ]]; then
    echo "Expected ${EXPECTED_SHARDS} non-empty ${table} files; found ${count}" >&2
    exit 2
  fi
done
archive_count="$(count_nonempty "${RUN_ROOT}/run-archives" '*.tar.zst')"
checksum_count="$(count_nonempty "${RUN_ROOT}/run-archives" '*.sha256')"
if [[ "${archive_count}" -ne "${EXPECTED_SHARDS}" || "${checksum_count}" -ne "${EXPECTED_SHARDS}" ]]; then
  echo "Archive/checksum count mismatch: ${archive_count}/${checksum_count}" >&2
  exit 2
fi
if find "${RUN_ROOT}" -type d -name runs -print -quit | grep -q .; then
  echo "Durable runs directory remains under ${RUN_ROOT}" >&2
  exit 2
fi
if find "${INPUT_LOG_DIR}" -type f -name "fdr_${ARRAY_JOB_ID}_*.err" -size +0c -print -quit | grep -q .; then
  echo "Non-empty calibration stderr exists for array ${ARRAY_JOB_ID}" >&2
  exit 2
fi

(
  cd "${RUN_ROOT}/run-archives"
  for checksum in ./*.sha256; do
    sha256sum -c "${checksum}" >/dev/null
  done
)
mapfile -t representative_archives < <(
  find "${RUN_ROOT}/run-archives" -maxdepth 1 -type f -name '*.tar.zst' | sort | awk \
    -v expected="${EXPECTED_SHARDS}" \
    'NR==1 || (expected > 2 && NR==int((expected+1)/2)) || NR==expected'
)
representative_count=3
if [[ "${EXPECTED_SHARDS}" -lt "${representative_count}" ]]; then
  representative_count="${EXPECTED_SHARDS}"
fi
if [[ "${#representative_archives[@]}" -ne "${representative_count}" ]]; then
  echo "Could not select ${representative_count} representative archives" >&2
  exit 2
fi
for archive in "${representative_archives[@]}"; do
  zstd -q -t "${archive}"
  tar -I zstd -tf "${archive}" >/dev/null
done

PYTHONPATH="${PROJECT_DIR}/src${PYTHONPATH:+:${PYTHONPATH}}" "${PYTHON_BIN}" \
  "${PROJECT_DIR}/experiments/002_relative_weighting/fdr_calibration_postprocess.py" \
  --root "${RUN_ROOT}" \
  --outdir "${OUTPUT_ROOT}" \
  --expected-shards "${EXPECTED_SHARDS}" \
  --expected-reps "${EXPECTED_REPS}" \
  --training-reps "${TRAINING_REPS}"

truth_tsv="${OUTPUT_ROOT}/truth_branch_audit.tsv"
truth_archive="${truth_tsv}.zst"
zstd -q -T"${ZSTD_THREADS}" -"${ZSTD_LEVEL}" "${truth_tsv}" -o "${truth_archive}"
zstd -q -t "${truth_archive}"
rm -f "${truth_tsv}"

manifest_tmp="${OUTPUT_ROOT}/.postprocess_manifest.tmp.${SLURM_JOB_ID:-$$}"
{
  printf 'field\tvalue\n'
  printf 'completed_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'postprocess_source_commit\t%s\n' "${SOURCE_ID}"
  printf 'input_array_job\t%s\n' "${ARRAY_JOB_ID}"
  printf 'input_run_root\t%s\n' "${RUN_ROOT}"
  printf 'input_manifest_sha256\t%s\n' "$(sha256sum "${RUN_ROOT}/run_manifest.tsv" | awk '{print $1}')"
  printf 'expected_shards\t%s\n' "${EXPECTED_SHARDS}"
  printf 'expected_replicates\t%s\n' "${EXPECTED_REPS}"
  printf 'training_replicates\t%s\n' "${TRAINING_REPS}"
  printf 'truth_archive_sha256\t%s\n' "$(sha256sum "${truth_archive}" | awk '{print $1}')"
} > "${manifest_tmp}"
mv "${manifest_tmp}" "${OUTPUT_ROOT}/postprocess_manifest.tsv"
(
  cd "${OUTPUT_ROOT}"
  find . -maxdepth 1 -type f ! -name postprocess_sha256.txt -print0 | sort -z | \
    xargs -0 sha256sum > postprocess_sha256.txt
  sha256sum -c postprocess_sha256.txt >/dev/null
)
if find "${OUTPUT_ROOT}" -maxdepth 1 -type f -name '*.tmp.*' -print -quit | grep -q .; then
  echo "Temporary files remain after postprocessing" >&2
  exit 2
fi
printf 'postprocess_output=%s\n' "${OUTPUT_ROOT}"
printf 'truth_archive=%s\n' "${truth_archive}"
printf 'checksums=verified\n'
