# Confirmation of the formulation experiment (TM-M2StructAlign)

Does optimising the lDDT alone (a genetic algorithm) give better alignments than the published
formulation (NSGA-II with SOPwTP and lDDT)? A preliminary experiment on three development datasets
said yes in two of them; this campaign confirms it on datasets not used before and checks a larger
budget. The protocol was fixed before running it.

| Cells | Dataset | Evaluations | Purpose |
| --- | --- | --- | --- |
| 1–3 | classB1, classC_19, classA_002_19 | 10,000 | Confirmation on unused datasets |
| 4–6 | classF, classT2_19, classA_001_19 | 20,000 | Does the advantage hold with a larger budget? |

| | |
| --- | --- |
| Configurations | `ga` (`--algorithm GENETIC_ALGORITHM --objectives LDDT`) and `nsgaii` (defaults) |
| Runs | seeds 1–6 × 2 configurations = 12 per cell, 72 in all |
| Common | population 100, seed alignments of clustalw, kalign, mafft, muscle5, tcoffee and tmaligner, 4 evaluation threads per run |
| Primary variable | Pair F1 against GPCRdb of the alignment with the best lDDT (rule R1) |
| Hypothesis | the GA beats NSGA-II in at least 2 of cells 1–3 with 5 or 6 seeds in favour (exact paired Wilcoxon) |
| Resources per task | 50 cores (12 runs at a time), 64 GB, 4 h limit |

## Contents

| Path | What |
| --- | --- |
| `TM-MSAStructAligner-1.0-jar-with-dependencies.jar` | The tool, with jMetal and Evolver inside |
| `TOOL_COMMIT`, `JMETAL_COMMIT`, `EVOLVER_COMMIT`, `DATASET_VERSION` | Provenance of the jar and the data |
| `source/` | The source code of the tool at `TOOL_COMMIT` |
| `data/GPCRdb/` | Only what the runs and the evaluation read, for the six datasets |
| `cells.tsv`, `configurations.tsv`, `campaign.env` | The campaign: cells, configurations and parameters |
| `slurm_campaign.sh` | One cell per array task |
| `submit.sh`, `check.sh`, `analyse.sh` | Submit, follow and analyse |
| `scripts/` | `formulation_study.py` and `compare_runs.py` (used by `analyse.sh`) |

## Steps

1. Copy this directory (or unzip the archive) to the cluster.
2. `./submit.sh` submits the six cells as one job array, then returns: the jobs run in the
   background. `./submit.sh 3` or `./submit.sh 1,4-6` submits some cells only.
3. Follow them with `squeue -u $USER` and `./check.sh`. `check.sh` can be run at any time; it
   lists the runs that are finished, failed or pending in each cell.
4. If a cell did not finish (time limit, node failure, a failed run), run `./check.sh` once its job
   has ended and the command it prints, e.g. `./submit.sh 2`. A resubmitted cell only runs again
   the runs without `metadata.json`; the finished ones are never repeated.
5. Copy `results/` (and `logs/` if something failed) back, and run `./analyse.sh` where Java and
   Python 3 are available. It writes `analysis/<cell>.csv` and `analysis/<cell>_comparison.txt`.

## Results

```text
results/<cell>/ejecuciones/<dataset>/<configuration>_seed_<n>/   one run (FUN.csv, MSASol*.fasta,
                                                                  config.json, metadata.json, ...)
results/<cell>/failed/<run>/      the files of a run that failed, if any
results/<cell>/status.tsv         run, exit code, time and job of every attempt
results/<cell>/progress.txt       start and end of every run
logs/<job>_<array>_<task>.out     log of a task (node, cores, versions, summary)
logs/runs/<cell>/<run>.log        output of every run
```

Each run writes its configuration and provenance (`config.json`, `metadata.json`: commit of the
tool, SHA-256 of the seed alignments, mapping of the structures, complete configuration of the
algorithm), so it can be repeated exactly with the same jar, data and seed: the number of threads
does not change the results.

## Estimated cost

Measured on a laptop with 4 evaluation threads per run: classF takes about 40 s with 10,000
evaluations; classT2_19 and classA_001_19 about 70 s. classC_19 (long proteins) is expected to be
the slowest (a few minutes). A cell takes as long as its slowest run, so each task should take
between 2 and 15 minutes, and the whole campaign about 10–20 core-hours. The 4 h limit leaves a
wide margin.
