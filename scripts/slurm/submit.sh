#!/usr/bin/env bash
#
# Submits the cells of the campaign as one slurm job array (one task per line of cells.tsv).
#
#   ./submit.sh                 every cell
#   ./submit.sh 2,5             the given cells (the --array specification of slurm)
#
# It only submits and returns: the jobs run in the background. Follow them with squeue and
# ./check.sh. It creates logs/ first, because sbatch redirects the output of a job there before
# the job starts.

set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"
source campaign.env

if [ ! -f "$JAR" ] || [ ! -d data/GPCRdb ]; then
  echo "The jar or data/GPCRdb are missing: build the campaign with prepare_campaign.sh" >&2
  exit 1
fi

CELLS=$(awk -F'\t' 'NR > 1 && $0 !~ /^#/' cells.tsv | wc -l | tr -d ' ')
ARRAY=${1:-1-$CELLS}
mkdir -p logs
echo "Submitting cells $ARRAY of $CELLS ($CAMPAIGN)"
sbatch --job-name="tm_$CAMPAIGN" --array="$ARRAY" slurm_campaign.sh
