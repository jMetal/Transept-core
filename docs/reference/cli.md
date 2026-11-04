# Command-line reference

The JAR built by `mvn clean package` runs the command line:

```bash
java -jar target/TM-MSAStructAligner-1.0-jar-with-dependencies.jar <command> [options]
```

| Command | Purpose |
| --- | --- |
| `run` | Align one dataset and write the final population |
| `baseline` | Score the precomputed alignments of one or more datasets |

`--help` prints the options of a command and `--version` the version. Options can be stored in a
file, one per line, and passed as `@file.args`.

## `run`

```bash
java -jar target/TM-MSAStructAligner-1.0-jar-with-dependencies.jar run \
  --data-root ../Dataset_Structural_GPCRs/GPCRdb \
  --dataset classF \
  --seeds clustalw,kalign,mafft,muscle5,tcoffee,tmaligner \
  --output-root results --run-id seed_000 --seed 0
```

### Input and output

| Option | Default | Meaning |
| --- | --- | --- |
| `--data-root` | (required) | `GPCRdb` directory of the [dataset repository](https://github.com/TM-M2StructAlign/Dataset_Structural_GPCRs) |
| `--dataset` | (required) | Dataset name, e.g. `classF` or `classA_001_19` |
| `--seeds` | (required) | Tools whose `precomputed/<dataset>/<tool>_<dataset>.fasta` seed the initial population, comma-separated |
| `--output-root` | (required) | Root of the results tree |
| `--run-id` | (required) | Run identifier; the run writes to `<output-root>/ejecuciones/<dataset>/<run-id>/`, which must not exist |
| `--seed` | (required) | Seed of the random generator; the same seed and options reproduce the run exactly |
| `--progress-interval` | 100 | Evaluations between two lines of `progress.jsonl` |

### Search

| Option | Default | Meaning |
| --- | --- | --- |
| `--algorithm` | `NSGA_II` | Search algorithm (see below) |
| `--max-evaluations` | 10000 | Stopping criterion |
| `--population-size` | 100 | Population size; also the offspring size and, for PAES, the archive size |
| `--threads` | 1 | Parallel evaluation threads; results do not depend on it |
| `--crossover-probability` | 0.8 | Probability of applying the crossover (single point, never cuts a transmembrane helix) to each child |
| `--mutation` | `TM_SAFE_GAP_SHIFT` | Mutation operator (see below) |
| `--mutation-probability` | 0.2 (1.0 for PAES) | Probability of mutating each child |

The initial population holds the distinct seeds and, for the rest, distinct crossover children
([T2, Step 7](../tutorials/T2-problem-formulation.md#step-7-the-initial-population)).

| `--algorithm` | Algorithm | Notes |
| --- | --- | --- |
| `NSGA_II` | NSGA-II | Published configuration |
| `NSGA_III` | NSGA-III | The population size must be a multiple of 4 |
| `AGE_MOEA`, `AGE_MOEA_II` | AGE-MOEA, AGE-MOEA-II | |
| `SMS_EMOA` | SMS-EMOA | Steady state: `--threads` does not speed it up |
| `MOEA_D` | MOEA/D | Normalised PBI aggregation, neighbourhood of 20, seeds spread over the subproblems; steady state |
| `PAES` | PAES, (1+1) with a crowding-distance archive | Mutation only; starts from the first seed |
| `RANDOM_SEARCH` | Random search, gaps at random | Sanity check (weak baseline); ignores the seeds and the operators |
| `RANDOM_SEEDED` | Random search, seeds + 10 random mutations | Sanity check (strong baseline) |

| `--mutation` | Never breaks a helix | Effect |
| --- | --- | --- |
| `TM_SAFE_GAP_SHIFT` (default) | Yes | Moves a gap run of one row to a position that does not split a helix |
| `MERGE_GAP_RUNS` | Yes | Joins two consecutive gap runs of one row |
| `RANDOM_CHOICE` | Yes | `TM_SAFE_GAP_SHIFT` or `MERGE_GAP_RUNS`, chosen uniformly each time |

The [formulation tutorial](../tutorials/T2-problem-formulation.md) describes every operator.

### Objectives and structures

| Option | Default | Meaning |
| --- | --- | --- |
| `--mapping-policy` | `VERIFIED` | How AlphaFold residues are mapped to the sequence: `VERIFIED` (numbering checked by residue identity, sequence alignment otherwise) or `NUMBERED_ONLY` (numbering only, structure dropped if residues disagree) |

## `baseline`

Scores every precomputed alignment of the given datasets with the two objectives.

```bash
java -jar target/TM-MSAStructAligner-1.0-jar-with-dependencies.jar baseline \
  --data-root ../Dataset_Structural_GPCRs/GPCRdb \
  --datasets classF,classT2_19 \
  --output baselines.csv
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--data-root` | (required) | `GPCRdb` directory |
| `--datasets` | (required) | Dataset names, comma-separated |
| `--output` | standard output | CSV file to write (see [outputs](outputs.md)) |
| `--mapping-policy` | as in `run` | |

## Comparison with the reference alignments

`EvaluationMain` compares every alignment of a results tree with a GPCRdb reference alignment
(see [metrics](metrics.md)). It is not a subcommand of the JAR's command line yet:

```bash
java -cp target/TM-MSAStructAligner-1.0-jar-with-dependencies.jar \
  org.tmm2structalign.evaluation.EvaluationMain \
  --reference ../Dataset_Structural_GPCRs/GPCRdb/reference_alignments/classF.fasta \
  --topology ../Dataset_Structural_GPCRs/GPCRdb/sequences/tmregions/classF_predicted_topologies.3line \
  --msa-root results/ejecuciones/classF \
  --output results/evaluation_classF.csv
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--reference` | (required) | Reference alignment (FASTA) |
| `--topology` | (required) | Three-line topology file of the dataset |
| `--msa-root` | (required) | Directory searched recursively for alignments |
| `--glob` | `MSA*.fasta` | Pattern of the alignment files |
| `--output` | `<msa-root>/evaluation_metrics.csv` | CSV report; it never overwrites an input file |
