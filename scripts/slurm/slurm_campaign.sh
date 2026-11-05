#!/usr/bin/env bash
#
# One cell of a TM-M2StructAlign campaign per array task: every configuration of
# configurations.tsv with every seed of campaign.env, on the dataset and budget of the cell
# (line <task id> of cells.tsv). The runs of a task are launched at the same time
# (PARALLEL_RUNS runs of THREADS_PER_RUN threads each) on the local scratch of the node.
#
# Meant to run unattended, in the background:
#   - a failed run does not stop the others; its exit code goes to results/<cell>/status.tsv and
#     its files to results/<cell>/failed/<run-id>/, and the task ends with an error;
#   - every finished run is copied back as soon as it ends, so a task killed by the time limit or
#     the node keeps what it finished;
#   - a resubmitted cell only launches the runs without metadata.json (see check.sh);
#   - logs/runs/<cell>/<run-id>.log has the output of every run, results/<cell>/progress.txt the
#     start and end of every run;
#   - the local scratch is always removed, also on failure.
#
# Submit it with submit.sh, from the directory of the campaign (prepare_campaign.sh builds it).

#SBATCH -J tm_campaign
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=64gb
#SBATCH --time=04:00:00
#SBATCH --constraint=cal
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=you@example.org

set -uo pipefail

HERE=${SLURM_SUBMIT_DIR:-$PWD}
cd "$HERE" || exit 1
# shellcheck source=/dev/null
source campaign.env

if command -v module > /dev/null 2>&1; then
  module load "${JAVA_MODULE:-java/jdk-25}"
fi

TASK=${SLURM_ARRAY_TASK_ID:?Run it as a slurm array task (or set SLURM_ARRAY_TASK_ID)}
LINE=$(awk -F'\t' -v n="$TASK" 'NR > 1 && $0 !~ /^#/ && ++i == n' cells.tsv)
if [ -z "$LINE" ]; then
  echo "cells.tsv has no cell $TASK" >&2
  exit 1
fi
IFS=$'\t' read -r CELL DATASET EVALUATIONS <<< "$LINE"
EVALUATIONS=${EVALUATIONS_OVERRIDE:-$EVALUATIONS}

SCRATCH="${LOCALSCRATCH:-${TMPDIR:-/tmp}}/$USER/${SLURM_JOB_ID:-local}_$TASK"
trap 'rm -rf --one-file-system "${SCRATCH:?}" 2> /dev/null || rm -rf "${SCRATCH:?}"' EXIT

RESULTS="$HERE/results/$CELL"
mkdir -p "$RESULTS/ejecuciones/$DATASET" "$RESULTS/failed" "$HERE/logs/runs/$CELL"
touch "$RESULTS/status.tsv"

echo "cell $CELL ($DATASET, $EVALUATIONS evaluations) on $(hostname), $(nproc 2> /dev/null || getconf _NPROCESSORS_ONLN) cores, $(date)"
echo "tool $(cat TOOL_COMMIT 2> /dev/null), jMetal $(cat JMETAL_COMMIT 2> /dev/null)," \
  "Evolver $(cat EVOLVER_COMMIT 2> /dev/null), dataset $(cat DATASET_VERSION 2> /dev/null)"

# Copy the jar and only the files of this dataset to the local scratch
DATA="$SCRATCH/GPCRdb"
mkdir -p "$DATA/sequences/tmregions" "$DATA/distances" "$DATA/alphafold_pdb" "$DATA/precomputed"
cp "$JAR" "$SCRATCH/" || exit 1
cp "data/GPCRdb/sequences/tmregions/${DATASET}_predicted_topologies.3line" \
  "$DATA/sequences/tmregions/" || exit 1
cp -R "data/GPCRdb/distances/$DATASET" "$DATA/distances/" || exit 1
cp -R "data/GPCRdb/precomputed/$DATASET" "$DATA/precomputed/" || exit 1
if [ -d "data/GPCRdb/alphafold_pdb/$DATASET" ]; then
  cp -R "data/GPCRdb/alphafold_pdb/$DATASET" "$DATA/alphafold_pdb/"
fi

# The runs of this cell that have not finished yet: "<configuration> <seed>" per line
PENDING="$SCRATCH/pending.txt"
: > "$PENDING"
while IFS=$'\t' read -r NAME _; do
  for SEED in $SEEDS; do
    if [ ! -f "$RESULTS/ejecuciones/$DATASET/${NAME}_seed_${SEED}/metadata.json" ]; then
      echo "$NAME $SEED" >> "$PENDING"
    fi
  done
done < <(awk -F'\t' 'NR > 1 && $0 !~ /^#/' configurations.tsv)
echo "$(wc -l < "$PENDING" | tr -d ' ') runs to do"

run_one() {
  local name=$1 seed=$2 options run_id run_dir log code
  options=$(awk -F'\t' -v n="$name" 'NR > 1 && $1 == n { print $2 }' "$HERE/configurations.tsv")
  run_id="${name}_seed_${seed}"
  log="$HERE/logs/runs/$CELL/$run_id.log"
  echo "$(date '+%F %T') start $run_id" >> "$RESULTS/progress.txt"
  # shellcheck disable=SC2086
  java -Xmx"$JAVA_HEAP" -jar "$SCRATCH/$(basename "$JAR")" run \
    --data-root "$DATA" --dataset "$DATASET" --seeds "$SEED_METHODS" \
    --output-root "$SCRATCH/out" --run-id "$run_id" --seed "$seed" \
    --max-evaluations "$EVALUATIONS" --population-size "$POPULATION" \
    --threads "$THREADS_PER_RUN" $options > "$log" 2>&1
  code=$?
  run_dir="$SCRATCH/out/ejecuciones/$DATASET/$run_id"
  if [ "$code" -eq 0 ] && [ -f "$run_dir/metadata.json" ]; then
    rm -rf "${RESULTS:?}/ejecuciones/$DATASET/$run_id" "${RESULTS:?}/failed/$run_id"
    cp -R "$run_dir" "$RESULTS/ejecuciones/$DATASET/"
  else
    code=${code/#0/1}
    rm -rf "${RESULTS:?}/failed/$run_id"
    if [ -d "$run_dir" ]; then cp -R "$run_dir" "$RESULTS/failed/"; fi
  fi
  printf '%s\t%s\t%s\t%s\n' "$run_id" "$code" "$(date '+%F %T')" "${SLURM_JOB_ID:-local}" \
    >> "$RESULTS/status.tsv"
  echo "$(date '+%F %T') end $run_id (exit $code)" >> "$RESULTS/progress.txt"
}
export -f run_one
export HERE CELL DATASET EVALUATIONS SCRATCH RESULTS DATA JAR JAVA_HEAP SEED_METHODS POPULATION \
  THREADS_PER_RUN

tr '\n' '\0' < "$PENDING" \
  | xargs -0 -P "$PARALLEL_RUNS" -I{} bash -c 'run_one $1' _ {}

# Exit with an error if a run of this cell is still not finished
FAILED=0
while IFS=$'\t' read -r NAME _; do
  for SEED in $SEEDS; do
    if [ ! -f "$RESULTS/ejecuciones/$DATASET/${NAME}_seed_${SEED}/metadata.json" ]; then
      echo "not finished: ${NAME}_seed_${SEED} (see logs/runs/$CELL/${NAME}_seed_${SEED}.log)"
      FAILED=$((FAILED + 1))
    fi
  done
done < <(awk -F'\t' 'NR > 1 && $0 !~ /^#/' configurations.tsv)
echo "cell $CELL done, $FAILED runs not finished, $(date)"
[ "$FAILED" -eq 0 ]
