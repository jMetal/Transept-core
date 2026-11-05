#!/usr/bin/env bash
#
# Evaluates every run against the GPCRdb reference and compares the configurations, cell by
# cell, with the indicators of ANALYSIS_INDICATORS (campaign.env) against the configuration
# ANALYSIS_BASELINE. It needs Java and Python 3 (standard library only). Run it where results/
# is, usually after copying it back from the cluster.
#
#   ./analyse.sh

set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"
source campaign.env

mkdir -p analysis
while IFS=$'\t' read -r CELL DATASET EVALUATIONS; do
  if [ ! -d "results/$CELL/ejecuciones/$DATASET" ]; then
    echo "results/$CELL has no runs yet" >&2
    continue
  fi
  python3 scripts/formulation_study.py "results/$CELL" "$DATASET" --data-root data/GPCRdb \
    --output "analysis/$CELL.csv" --jar "$JAR"
  {
    echo "# $CELL: $DATASET, $EVALUATIONS evaluations"
    for INDICATOR in $ANALYSIS_INDICATORS; do
      python3 scripts/compare_runs.py "analysis/$CELL.csv" --baseline "$ANALYSIS_BASELINE" \
        --indicator "$INDICATOR"
      echo
    done
  } | tee "analysis/${CELL}_comparison.txt"
done < <(awk -F'\t' 'NR > 1 && $0 !~ /^#/' cells.tsv)
