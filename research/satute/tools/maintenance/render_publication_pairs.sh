#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: render_publication_pairs.sh --summary FILE --detail FILE --out-root DIR
       [--suite gtr|protein|misspec|all] [--expected-reps N]
       [--bootstrap-reps N] [--render-script FILE] [--force]

Each pair/rule directory must contain the R renderer's 16 non-empty outputs:
three 2D PDFs/SVGs/PNGs, three vector 3D PDFs, three self-contained HTML
files, and one paired-difference TSV.
EOF
}

summary=""
detail=""
out_root=""
suite="all"
expected_reps=1000
bootstrap_reps=5000
render_script=""
force=0

while (($#)); do
  case "$1" in
    --summary) summary=$2; shift 2 ;;
    --detail) detail=$2; shift 2 ;;
    --out-root) out_root=$2; shift 2 ;;
    --suite) suite=$2; shift 2 ;;
    --expected-reps) expected_reps=$2; shift 2 ;;
    --bootstrap-reps) bootstrap_reps=$2; shift 2 ;;
    --render-script) render_script=$2; shift 2 ;;
    --force) force=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z "$summary" || -z "$detail" || -z "$out_root" ]]; then
  usage >&2
  exit 2
fi
[[ -s "$summary" ]] || { echo "Missing or empty summary: $summary" >&2; exit 1; }
[[ -s "$detail" ]] || { echo "Missing or empty detail: $detail" >&2; exit 1; }
if [[ -z "$render_script" ]]; then
  render_script="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../experiments/002_relative_weighting" && pwd)/render_formula_comparison.R"
fi
[[ -s "$render_script" ]] || { echo "Missing renderer: $render_script" >&2; exit 1; }

case "$suite" in
  gtr)
    pairs=(
      GTR_PF06346:GTR_PF06346
      GTR_PF06346_G4:GTR_PF06346_G4
      GTR_PF06346_I_G4:GTR_PF06346_I_G4
    )
    ;;
  protein)
    pairs=(LG_G4:LG_G4 WAG_G4:WAG_G4 JTT_G4:JTT_G4 Q.PFAM_G4:Q.PFAM_G4)
    ;;
  misspec)
    pairs=(GTR_PF06346:JC GTR_PF06346:K2P GTR_PF06346:F81)
    ;;
  all)
    pairs=(
      GTR_PF06346:GTR_PF06346
      GTR_PF06346_G4:GTR_PF06346_G4
      GTR_PF06346_I_G4:GTR_PF06346_I_G4
      LG_G4:LG_G4 WAG_G4:WAG_G4 JTT_G4:JTT_G4 Q.PFAM_G4:Q.PFAM_G4
      GTR_PF06346:JC GTR_PF06346:K2P GTR_PF06346:F81
    )
    ;;
  *) echo "Invalid suite: $suite" >&2; exit 2 ;;
esac

rules=(unadjusted taxon_bonferroni by_fdr)
for pair in "${pairs[@]}"; do
  simulation_model=${pair%%:*}
  evaluation_model=${pair##*:}
  for rule in "${rules[@]}"; do
    output_dir="$out_root/$rule/${simulation_model}__${evaluation_model}"
    marker="$output_dir/paired_formula_difference_rule_${rule}__sim_${simulation_model}__eval_${evaluation_model}.tsv"
    if ((force == 0)) && [[ -s "$marker" ]]; then
      echo "SKIP $simulation_model:$evaluation_model $rule (marker exists)"
      continue
    fi
    mkdir -p "$output_dir"
    echo "RENDER $simulation_model:$evaluation_model $rule"
    Rscript "$render_script" \
      --summary "$summary" --detail "$detail" --outdir "$output_dir" \
      --simulation-model "$simulation_model" --evaluation-model "$evaluation_model" \
      --expected-reps "$expected_reps" --bootstrap-reps "$bootstrap_reps" \
      --decision-rule "$rule"
    count=$(find "$output_dir" -maxdepth 1 -type f -size +0c | wc -l | tr -d ' ')
    [[ "$count" == 16 ]] || {
      echo "Output contract failed for $output_dir: expected 16 non-empty files, got $count" >&2
      exit 1
    }
  done
done

echo "Publication pair rendering passed: suite=$suite out_root=$out_root"
