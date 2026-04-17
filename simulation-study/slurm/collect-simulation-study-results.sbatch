#!/usr/bin/env bash
#SBATCH --partition=standard
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

set -euo pipefail

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <repo-root> <output-dir>" >&2
  exit 1
fi

REPO_ROOT="$1"
OUTPUT_DIR="$2"

cd "$REPO_ROOT"
Rscript simulation-study/scripts/collect-simulation-study-results.R \
  --output-dir="$OUTPUT_DIR"
