#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
ROOT_DIR="${ROOT_DIR:-$(cd "${PROJECT_DIR}/../.." && pwd)}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
export PYTHONPATH="${PROJECT_DIR}/src${PYTHONPATH:+:${PYTHONPATH}}"
SHARD_ROOT="${SHARD_ROOT:-/work/satute-head-to-head}"
OUTDIR="${OUTDIR:-${PROJECT_DIR}/artifacts/cluster/aws/full-run}"
FIGURE_DIR="${FIGURE_DIR:-${OUTDIR}/figures}"
EXPECTED_REPS="${EXPECTED_REPS:-1000}"
TREE_CASES="${TREE_CASES:-five_external,sixteen_internal}"
PLOT_MODEL_PAIRS="${PLOT_MODEL_PAIRS:-JC:JC}"
PLOT_SCENARIO_SET="${PLOT_SCENARIO_SET:-fig2}"
PLOT_DECISION_RULES="${PLOT_DECISION_RULES:-unadjusted,taxon_bonferroni,by_fdr}"
VERIFY_Z_TOLERANCE="${VERIFY_Z_TOLERANCE:-1e-4}"
ALLOW_INCOMPLETE="${ALLOW_INCOMPLETE:-0}"

render_extra_args=()
case "${ALLOW_INCOMPLETE}" in
  0|false) ;;
  1|true) render_extra_args=(--allow-incomplete) ;;
  *)
    echo "ALLOW_INCOMPLETE must be 0, 1, false, or true; got: ${ALLOW_INCOMPLETE}" >&2
    exit 2
    ;;
esac

DETAILS=()
while IFS= read -r detail_path; do
  DETAILS+=("${detail_path}")
done < <(find "${SHARD_ROOT}" -path '*/head_to_head_detail.tsv' -type f | sort)
if [[ "${#DETAILS[@]}" -eq 0 ]]; then
  echo "No shard detail TSV files found under ${SHARD_ROOT}" >&2
  exit 1
fi

mkdir -p "${OUTDIR}" "${FIGURE_DIR}"

"${PYTHON_BIN}" "${PROJECT_DIR}/experiments/002_relative_weighting/combine.py" \
  --inputs "${DETAILS[@]}" \
  --outdir "${OUTDIR}"

"${PYTHON_BIN}" "${PROJECT_DIR}/tests/unit/test_head_to_head_calculations.py" \
  --detail "${OUTDIR}/head_to_head_detail.tsv" \
  --z-tolerance "${VERIFY_Z_TOLERANCE}"

"${PYTHON_BIN}" "${PROJECT_DIR}/experiments/002_relative_weighting/profile.py" \
  --timing "${OUTDIR}/head_to_head_timing.tsv"

IFS=',' read -r -a plot_model_pairs <<< "${PLOT_MODEL_PAIRS}"
for pair in "${plot_model_pairs[@]}"; do
  simulation_model="${pair%%:*}"
  evaluation_model="${pair#*:}"
  if [[ -z "${simulation_model}" || -z "${evaluation_model}" || "${simulation_model}" == "${pair}" ]]; then
    echo "Invalid PLOT_MODEL_PAIRS item '${pair}'. Expected SIMULATION_MODEL:EVALUATION_MODEL." >&2
    exit 2
  fi
  IFS=',' read -r -a decision_rules <<< "${PLOT_DECISION_RULES}"
  for decision_rule in "${decision_rules[@]}"; do
    "${PYTHON_BIN}" "${PROJECT_DIR}/experiments/002_relative_weighting/render.py" \
      --summary "${OUTDIR}/head_to_head_summary.tsv" \
      --outdir "${FIGURE_DIR}" \
      --expected-reps "${EXPECTED_REPS}" \
      --tree-cases "${TREE_CASES}" \
      --simulation-model "${simulation_model}" \
      --evaluation-model "${evaluation_model}" \
      --scenario-set "${PLOT_SCENARIO_SET}" \
      --decision-rule "${decision_rule}" \
      "${render_extra_args[@]}"
  done
done

echo "Combined results: ${OUTDIR}"
echo "Full-run figures: ${FIGURE_DIR}"
