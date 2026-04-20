#!/usr/bin/env bash
#SBATCH --partition=standard
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <repo-root> <output-dir> <manifest-path>" >&2
  exit 1
fi

REPO_ROOT="$1"
OUTPUT_DIR="$2"
MANIFEST_PATH="$3"
TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

cd "$REPO_ROOT"
Rscript simulation-study/scripts/run-ci-study-cell.R \
  --output-dir="$OUTPUT_DIR" \
  --manifest="$MANIFEST_PATH" \
  --task-id="$TASK_ID"
