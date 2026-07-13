#!/usr/bin/env bash
#SBATCH --job-name=satute-full-h2h
#SBATCH --account=menche
#SBATCH --array=0-999%64
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=04:00:00
#SBATCH --output=research/satute/artifacts/cluster/lisc/logs/full_%A_%a.out
#SBATCH --error=research/satute/artifacts/cluster/lisc/logs/full_%A_%a.err

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "${SLURM_SUBMIT_DIR:-.}" && pwd)}"
PROJECT_DIR="${PROJECT_DIR:-${ROOT_DIR}/research/satute}"
WORK_DIR="${WORK_DIR:-/lisc/data/scratch/menche/sakalli/iq-tree-satute/full-head-to-head}"
SHARD_COUNT="${SHARD_COUNT:-1000}"

mkdir -p "${PROJECT_DIR}/artifacts/cluster/lisc/logs" "${WORK_DIR}"

ROOT_DIR="${ROOT_DIR}" PROJECT_DIR="${PROJECT_DIR}" WORK_DIR="${WORK_DIR}" SHARD_COUNT="${SHARD_COUNT}" \
  "${PROJECT_DIR}/tools/cluster/lisc/run_full_head_to_head_shard.sh"
