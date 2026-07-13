#!/usr/bin/env bash
#SBATCH --job-name=satute-protein-canary
#SBATCH --account=menche
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-/lisc/home/user/sakalli/projects/iq-tree-satute-protein-extension}"
PROJECT_DIR="${PROJECT_DIR:-${ROOT_DIR}/research/satute}"
WORK_DIR="${WORK_DIR:-/lisc/data/scratch/menche/sakalli/iq-tree-satute/full-head-to-head}"
IQTREE_BIN="${IQTREE_BIN:-${WORK_DIR}/iqtree3-build/iqtree3}"

export ROOT_DIR PROJECT_DIR WORK_DIR IQTREE_BIN
export RUN_NAME="${RUN_NAME:-protein-fixed-canary-${SLURM_JOB_ID}}"
export SIMULATOR=alisim
export PAPER_GRID=0
export SITE_LENGTHS=100
export BRANCH_LENGTHS=1.0
export REPS=1
export MODEL_PAIRS="LG:LG,WAG:WAG,JTT:JTT,Q.PFAM:Q.PFAM,LG_G4:LG_G4,WAG_G4:WAG_G4,JTT_G4:JTT_G4,Q.PFAM_G4:Q.PFAM_G4"
export SCENARIO_SET=misspecification
export TREE_CASES=five_external
export SHARD_INDEX=0
export SHARD_COUNT=1
export ARCHIVE_RUNS=0

"${PROJECT_DIR}/tools/cluster/lisc/run_full_head_to_head_shard.sh"
