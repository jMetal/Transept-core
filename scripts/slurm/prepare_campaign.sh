#!/usr/bin/env bash
#
# Builds the directory (and a zip of it) of a Transept campaign to run on a Slurm
# cluster, with only what it needs (clusters have disk quotas):
#
#   the jar with dependencies, and the commits of the tool, jMetal and Evolver it was built from
#   the source code of the tool at that commit (git archive)
#   data/GPCRdb/: for every dataset of cells.tsv, its topology, distances, AlphaFold models, the
#                 seed alignments of SEED_METHODS and its reference alignment
#   the campaign (cells.tsv, configurations.tsv, campaign.env, README.md) and the scripts
#
#   scripts/slurm/prepare_campaign.sh <campaign template> <output directory>
#   e.g. scripts/slurm/prepare_campaign.sh scripts/slurm/formulations ../campaigns/formulations
#
# Build the jar first (mvn package). The data, jMetal and Evolver are looked for next to this
# repository; override them with DATASET_REPOSITORY, JMETAL_REPOSITORY and EVOLVER_REPOSITORY.

set -euo pipefail

if [ $# -ne 2 ]; then
  echo "Usage: $0 <campaign template> <output directory>" >&2
  exit 2
fi
TEMPLATE=$(cd "$1" && pwd)
OUTPUT=$2
SCRIPTS=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
TOOL=$(cd "$SCRIPTS/../.." && pwd)
DATASET_REPOSITORY=${DATASET_REPOSITORY:-$TOOL/../Dataset_Structural_GPCRs}
JMETAL_REPOSITORY=${JMETAL_REPOSITORY:-$HOME/Softw/jMetal/jMetal}
EVOLVER_REPOSITORY=${EVOLVER_REPOSITORY:-$HOME/Softw/Evolver}
GPCRDB=$DATASET_REPOSITORY/GPCRdb

source "$TEMPLATE/campaign.env"

if [ ! -f "$TOOL/target/$JAR" ]; then
  echo "$TOOL/target/$JAR not found: run \"mvn package\" in $TOOL first" >&2
  exit 1
fi
if [ -e "$OUTPUT" ]; then
  echo "$OUTPUT exists: choose another directory or remove it" >&2
  exit 1
fi
for repository in "$TOOL" "$JMETAL_REPOSITORY" "$EVOLVER_REPOSITORY"; do
  if [ -n "$(git -C "$repository" status --porcelain -- src pom.xml '*/src' 2> /dev/null)" ]; then
    echo "Warning: $repository has uncommitted changes; its commit will not describe the jar" >&2
  fi
done

mkdir -p "$OUTPUT/source" "$OUTPUT/scripts"
cp "$TOOL/target/$JAR" "$OUTPUT/"
git -C "$TOOL" rev-parse HEAD > "$OUTPUT/TOOL_COMMIT"
git -C "$JMETAL_REPOSITORY" rev-parse HEAD > "$OUTPUT/JMETAL_COMMIT" 2> /dev/null \
  || echo unknown > "$OUTPUT/JMETAL_COMMIT"
git -C "$EVOLVER_REPOSITORY" rev-parse HEAD > "$OUTPUT/EVOLVER_COMMIT" 2> /dev/null \
  || echo unknown > "$OUTPUT/EVOLVER_COMMIT"
git -C "$DATASET_REPOSITORY" describe --tags --always --dirty > "$OUTPUT/DATASET_VERSION"
COMMIT=$(cut -c1-9 "$OUTPUT/TOOL_COMMIT")
git -C "$TOOL" archive --format=tar.gz --prefix="Transept-$COMMIT/" \
  -o "$OUTPUT/source/Transept-core-$COMMIT.tar.gz" HEAD

cp "$TEMPLATE"/cells.tsv "$TEMPLATE"/configurations.tsv "$TEMPLATE"/campaign.env \
  "$TEMPLATE"/README.md "$OUTPUT/"
cp "$SCRIPTS"/slurm_campaign.sh "$SCRIPTS"/submit.sh "$SCRIPTS"/check.sh "$SCRIPTS"/analyse.sh \
  "$OUTPUT/"
cp "$TOOL"/scripts/formulation_study.py "$TOOL"/scripts/compare_runs.py "$OUTPUT/scripts/"

DATA=$OUTPUT/data/GPCRdb
mkdir -p "$DATA/sequences/tmregions" "$DATA/distances" "$DATA/alphafold_pdb" \
  "$DATA/precomputed" "$DATA/reference_alignments"
for DATASET in $(awk -F'\t' 'NR > 1 && $0 !~ /^#/ { print $2 }' "$TEMPLATE/cells.tsv" | sort -u)
do
  cp "$GPCRDB/sequences/tmregions/${DATASET}_predicted_topologies.3line" \
    "$DATA/sequences/tmregions/"
  cp -R "$GPCRDB/distances/$DATASET" "$DATA/distances/"
  if [ -d "$GPCRDB/alphafold_pdb/$DATASET" ]; then
    cp -R "$GPCRDB/alphafold_pdb/$DATASET" "$DATA/alphafold_pdb/"
  fi
  mkdir -p "$DATA/precomputed/$DATASET"
  for METHOD in ${SEED_METHODS//,/ }; do
    cp "$GPCRDB/precomputed/$DATASET/${METHOD}_${DATASET}.fasta" "$DATA/precomputed/$DATASET/"
  done
  cp "$GPCRDB/reference_alignments/$DATASET.fasta" "$DATA/reference_alignments/"
done

(cd "$(dirname "$OUTPUT")" && zip -qr "$(basename "$OUTPUT").zip" "$(basename "$OUTPUT")")
echo "Ready: $OUTPUT ($(du -sh "$OUTPUT" | cut -f1)) and $OUTPUT.zip" \
  "($(du -sh "$OUTPUT.zip" | cut -f1)); tool $COMMIT, dataset $(cat "$OUTPUT/DATASET_VERSION")"
