#!/usr/bin/env bash
#SBATCH --job-name=satute-full-h2h
#SBATCH --account=menche
#SBATCH --array=0-767%64
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=72:00:00
#SBATCH --output=doc/satute-wiki/lisc/logs/full_%A_%a.out
#SBATCH --error=doc/satute-wiki/lisc/logs/full_%A_%a.err

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$(cd "${SLURM_SUBMIT_DIR:-.}" && pwd)}"
WORK_DIR="${WORK_DIR:-/lisc/data/scratch/menche/sakalli/iq-tree-satute/full-head-to-head}"
SHARD_COUNT="${SHARD_COUNT:-768}"

mkdir -p "${ROOT_DIR}/doc/satute-wiki/lisc/logs" "${WORK_DIR}"

ROOT_DIR="${ROOT_DIR}" WORK_DIR="${WORK_DIR}" SHARD_COUNT="${SHARD_COUNT}" \
  "${ROOT_DIR}/doc/satute-wiki/lisc/run_full_head_to_head_shard.sh"
