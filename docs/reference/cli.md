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
| `--max-evaluations` | 10000 | Stop after this many evaluations |
| `--max-time` | (none) | Stop after this many seconds instead; exclusive with `--max-evaluations`. Runs that stop by time are not reproducible exactly |
| `--population-size` | 100 | Population size; also the offspring size (unless the configuration sets it) and, for PAES, the archive size |
| `--threads` | 1 | Parallel evaluation threads; results do not depend on it |
| `--crossover-probability` | 0.8 | Probability of applying the crossover (single point, never cuts a transmembrane helix) to each child |
| `--mutation` | `TM_SAFE_GAP_SHIFT` | Mutation operator (see below) |
| `--mutation-probability` | 0.2 (1.0 for PAES) | Probability of mutating each child |
| `--configuration` | (none) | NSGA-II, SPEA2 and the GA: configuration from their parameter space, `--parameter value ...` ([configuration](configuration.md)). Unset parameters keep their default |
| `--configuration-file` | (none) | The same, read from the first non-blank line of a file |

`--mutation`, `--crossover-probability` and `--mutation-probability` set the corresponding
parameters of the configuration of the configurable algorithms; giving one of them and the same parameter in
`--configuration` is an error.

The initial population holds the distinct seeds and, for the rest, distinct crossover children
([T2, Step 7](../tutorials/T2-problem-formulation.md#step-7-the-initial-population)).

| `--algorithm` | Algorithm | Notes |
| --- | --- | --- |
| `NSGA_II` | NSGA-II | Published configuration; configurable ([configuration](configuration.md)) |
| `ASYNC_NSGA_II` | Asynchronous NSGA-II | Every evaluated solution enters the population at once, so all `--threads` stay busy; with more than one thread, runs do not repeat exactly |
| `SPEA2` | SPEA2-style (RDEMOEA with strength ranking and k-nearest-neighbour density) | The population plays the role of SPEA2's archive; configurable |
| `NSGA_III` | NSGA-III | The population size must be a multiple of 4 |
| `AGE_MOEA`, `AGE_MOEA_II` | AGE-MOEA, AGE-MOEA-II | |
| `SMS_EMOA` | SMS-EMOA | Steady state: `--threads` does not speed it up |
| `MOEA_D` | MOEA/D | Normalised PBI aggregation, neighbourhood of 20, seeds spread over the subproblems; steady state. With three objectives the population size must be one with weight vectors: 91, 100, 136, 300, 500, 600, 800, 1000 or 1200 |
| `PAES` | PAES, (1+1) with a crowding-distance archive | Mutation only; starts from the first seed |
| `RANDOM_SEARCH` | Random search, gaps at random | Sanity check (weak baseline); ignores the seeds and the operators |
| `RANDOM_SEEDED` | Random search, seeds + 10 random mutations | Sanity check (strong baseline) |
| `GENETIC_ALGORITHM` | Generational GA, binary tournament, (μ+λ) replacement | Single objective only (`--objectives LDDT`); same seeded population and operators as NSGA-II; configurable |
| `ASYNC_GENETIC_ALGORITHM` | Asynchronous steady-state GA | Single objective only; like `ASYNC_NSGA_II`, does not repeat exactly with more than one thread |

The multiobjective algorithms accept two or three objectives and the genetic algorithms one;
the random searches accept any number.

| `--mutation` | Never breaks a helix | Effect |
| --- | --- | --- |
| `TM_SAFE_GAP_SHIFT` (default) | Yes | Moves a gap run of one row to a position that does not split a helix |
| `MERGE_GAP_RUNS` | Yes | Joins two consecutive gap runs of one row |
| `RANDOM_CHOICE` | Yes | `TM_SAFE_GAP_SHIFT` or `MERGE_GAP_RUNS`, chosen uniformly each time |

The [formulation tutorial](../tutorials/T2-problem-formulation.md) describes every operator.

### Objectives and structures

| Option | Default | Meaning |
| --- | --- | --- |
| `--objectives` | `SOPWTP_LDDT` | Objectives to optimise, in the order of `FUN.csv`: `LDDT` (lDDT alone, for `GENETIC_ALGORITHM`), `SOPWTP_LDDT` (published formulation) or `SOPWTP_LDDT_COVERAGE` (adds the aligned-pair coverage, [metrics](metrics.md)) |
| `--mapping-policy` | `VERIFIED` | How AlphaFold residues are mapped to the sequence: `VERIFIED` (numbering checked by residue identity, sequence alignment otherwise) or `NUMBERED_ONLY` (numbering only, structure dropped if residues disagree) |

## `algorithms`

Prints, as JSON, every algorithm of `run --algorithm` with the numbers of objectives it accepts
(`minObjectives`, `maxObjectives`), its parameter space (`null` if it is not configurable) and
whether it is steady state, asynchronous and uses the crossover. The studio reads it instead of
keeping its own copy.

```bash
java -jar target/TM-MSAStructAligner-1.0-jar-with-dependencies.jar algorithms
```

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
