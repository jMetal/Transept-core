#!/usr/bin/env bash
#
# Reports, for every cell, the runs finished, failed and pending, and prints the command that
# resubmits the cells with unfinished runs. It can be used while the jobs run.
#
#   ./check.sh

set -uo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"
source campaign.env

TODO=()
TASK=0
while IFS=$'\t' read -r CELL DATASET EVALUATIONS; do
  TASK=$((TASK + 1))
  DONE=0 FAILED=0 PENDING=0
  MISSING=()
  while IFS=$'\t' read -r NAME _; do
    for SEED in $SEEDS; do
      RUN="${NAME}_seed_${SEED}"
      if [ -f "results/$CELL/ejecuciones/$DATASET/$RUN/metadata.json" ]; then
        DONE=$((DONE + 1))
      elif [ -d "results/$CELL/failed/$RUN" ] \
          || grep -q "^$RUN"$'\t' "results/$CELL/status.tsv" 2> /dev/null; then
        FAILED=$((FAILED + 1)); MISSING+=("$RUN (failed)")
      else
        PENDING=$((PENDING + 1)); MISSING+=("$RUN")
      fi
    done
  done < <(awk -F'\t' 'NR > 1 && $0 !~ /^#/' configurations.tsv)
  printf '%-4s %-16s %6s evaluations: %3d finished, %3d failed, %3d pending\n' \
    "$TASK" "$DATASET" "$EVALUATIONS" "$DONE" "$FAILED" "$PENDING"
  if [ ${#MISSING[@]} -gt 0 ]; then
    printf '       %s\n' "${MISSING[@]}" | head -12
    TODO+=("$TASK")
  fi
done < <(awk -F'\t' 'NR > 1 && $0 !~ /^#/' cells.tsv)

if [ ${#TODO[@]} -eq 0 ]; then
  echo "Every run has finished: ./analyse.sh"
else
  echo "Resubmit the unfinished cells (only their missing runs run again) once their jobs have"
  echo "ended (squeue): ./submit.sh $(IFS=,; echo "${TODO[*]}")"
fi
