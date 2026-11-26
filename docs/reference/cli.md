# Command-line reference

The JAR built by `mvn clean package` runs the command line:

```bash
java -jar target/transept-core-0.1-SNAPSHOT-jar-with-dependencies.jar <command> [options]
```

| Command | Purpose |
| --- | --- |
| `run` | Align one dataset and write the final population |
| `seeds` | Compute the internal seeds of datasets in advance |
| `algorithms` | Print the algorithms available, as JSON |
| `baseline` | Score the precomputed alignments of one or more datasets |
| `score` | Compare the alignments of a directory with a reference alignment |

`--help` prints the options of a command and `--version` the version. Options can be stored in a
file, one per line, and passed as `@file.args`.

## `run`

```bash
java -jar target/transept-core-0.1-SNAPSHOT-jar-with-dependencies.jar run \
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
| `--seeds` | `internal` | Seed alignments of the initial population, comma-separated: `internal` (up to 10 alignments computed by the tool with BioJava's progressive aligner, see the `seeds` command) and/or tools whose `precomputed/<dataset>/<tool>_<dataset>.fasta` exist |
| `--repair-seeds`, `--no-repair-seeds` | on | Move the gaps that split a transmembrane helix out of every seed before the search. The operators never create such gaps, so no alignment of the front has one. `--no-repair-seeds` keeps the gaps of the precomputed alignments; the internal seeds are always repaired |
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
| `--crossover` | `TM_SAFE_SINGLE_POINT` | Crossover operator (see below) |
| `--crossover-probability` | 0.8 | Probability of applying the crossover |
| `--cut-candidates` | 1 | Valid cuts the single-point crossover draws, of which it uses the one that adds the fewest gaps; 1 is a cut at random. Also applies to the crossover that fills the initial population, which is always the single-point one |
| `--mutation` | `TM_SAFE_GAP_SHIFT` | Mutation operator (see below) |
| `--mutation-probability` | 0.2 (1.0 for PAES) | Probability of mutating each child |
| `--max-shift` | 0 | Largest number of columns the mutation can move a gap run (`TM_SAFE_GAP_SHIFT` and `RANDOM_CHOICE`); 0 is no limit |
| `--configuration` | (none) | NSGA-II, SPEA2 and the GA: configuration from their parameter space, `--parameter value ...` ([configuration](configuration.md)). Unset parameters keep their default |
| `--configuration-file` | (none) | The same, read from the first non-blank line of a file |

#### Stopping a run early

Creating a file called `STOP` in the run directory (`<output-root>/ejecuciones/<dataset>/<run-id>/`)
ends the run at its next check, as if its budget were spent: it writes its results as usual, records
`"stoppedByUser": true` in `metadata.json` and deletes the file. The file is looked for at most
five times per second. If the algorithm has no results yet, it first finishes its first generation
(two populations evaluated), or its first solution in a random search. A run stopped this way cannot
be repeated exactly. The Monitor page of the Transept interface does this with its **Stop and keep
the results** button.

`--crossover`, `--mutation`, `--crossover-probability`, `--mutation-probability`, `--max-shift`
and `--cut-candidates` set the corresponding parameters of the configuration of the configurable algorithms; giving one of them and the same parameter in
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

| `--crossover` | Never cuts a helix | Effect |
| --- | --- | --- |
| `TM_SAFE_SINGLE_POINT` (default) | Yes | Takes the start of one parent and the rest of the other, padding with gaps to keep the second part aligned |
| `TM_SAFE_UNIFORM` | Yes | Exchanges blocks between the columns where both parents agree; adds no gap |

| `--mutation` | Never breaks a helix | Effect |
| --- | --- | --- |
| `TM_SAFE_GAP_SHIFT` (default) | Yes | Moves a gap run of one row to a position that does not split a helix |
| `MERGE_GAP_RUNS` | Yes | Joins two consecutive gap runs of one row |
| `RANDOM_CHOICE` | Yes | One of the other four mutations, chosen uniformly each time |
| `TM_SAFE_BLOCK_SHIFT` | Yes | Moves a block of gaps shared by several rows one column (by default) left or right |
| `TM_SAFE_WINDOW_REALIGN` | Yes | Realigns the residues of one row against the others inside a window of 50 columns (by default) |

The [formulation tutorial](../tutorials/T2-problem-formulation.md) describes every operator, and
the [operators reference](operators.md) lists their names and parameters in one place.

### Objectives and structures

| Option | Default | Meaning |
| --- | --- | --- |
| `--objectives` | `SOPWTP_LDDT` | Objectives to optimise, in the order of `FUN.csv`: `LDDT` (lDDT alone, for `GENETIC_ALGORITHM`), `SOPWTP_LDDT` (published formulation) or `SOPWTP_LDDT_COVERAGE` (adds the aligned-pair coverage, [metrics](metrics.md)) |
| `--mapping-policy` | `VERIFIED` | How AlphaFold residues are mapped to the sequence: `VERIFIED` (numbering checked by residue identity, sequence alignment otherwise) or `NUMBERED_ONLY` (numbering only, structure dropped if residues disagree) |

## `seeds`

Computes the internal seeds of datasets in advance, so that the runs with `--seeds internal` and
the same `--output-root` read them instead of computing them (a few seconds for 11–19 proteins,
minutes for the largest datasets).

```bash
java -jar target/transept-core-0.1-SNAPSHOT-jar-with-dependencies.jar seeds \
  --data-root ../Dataset_Structural_GPCRs/GPCRdb --datasets classF --output-root results
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--data-root`, `--datasets`, `--output-root` | (required) | As in `run`; the seeds go to `<output-root>/internal_seeds/<dataset>/` |
| `--count` | 10 | Distinct seeds per dataset, at most 10 |

The directory holds the alignments and `manifest.tsv`, with the fingerprint of the proteins
(identifiers, residues, topology) and the SHA-256 of every file. A run uses the stored seeds only
if they match its proteins and files; otherwise it computes them again and says why in its log.
The internal seeds are deterministic, so both ways give the same alignments. With very short or
almost identical sequences all the settings may give a single alignment; then add a seed of
another tool, because the population needs at least two distinct seeds.

## `algorithms`

Prints, as JSON, every algorithm of `run --algorithm` with the numbers of objectives it accepts
(`minObjectives`, `maxObjectives`), its parameter space (`null` if it is not configurable) and
whether it is steady state, asynchronous and uses the crossover.

```bash
java -jar target/transept-core-0.1-SNAPSHOT-jar-with-dependencies.jar algorithms
```

## `baseline`

Scores every precomputed alignment of the given datasets with the two objectives, or the alignment
files given with `--alignments`, such as the internal seeds that a run saved in its `seeds/`
folder (their paths are in its `metadata.json`).

```bash
java -jar target/transept-core-0.1-SNAPSHOT-jar-with-dependencies.jar baseline \
  --data-root ../Dataset_Structural_GPCRs/GPCRdb \
  --datasets classF,classT2_19 \
  --output baselines.csv
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--data-root` | (required) | `GPCRdb` directory |
| `--datasets` | (required) | Dataset names, comma-separated |
| `--alignments` | (none) | Aligned FASTA files to score instead of the precomputed ones, comma-separated; needs a single dataset. The method is the name of the file without `.fasta` |
| `--output` | standard output | CSV file to write (see [outputs](outputs.md)) |
| `--mapping-policy` | as in `run` | |

## `score`

Compares every alignment of a results tree with a reference alignment, such as those of GPCRdb
([metrics](metrics.md)). An alignment that cannot be compared does not stop the report: its row
says why ([outputs](outputs.md#evaluation-output)). The status is 0 if at least one alignment was
compared.

```bash
java -jar target/transept-core-0.1-SNAPSHOT-jar-with-dependencies.jar score \
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

## `repair`

Moves the gaps that split a transmembrane helix out of alignments, without any search: the same
repair that `run` applies to the seeds by default
([T8](../tutorials/T8-repairing-helix-gaps.md)). Each gap run inside a helix moves to the nearest
position outside the helices; the residues keep their order and the length of the alignment does
not change. The repaired alignments are written with the same file names, and the command prints how
many gap runs it moved in each one.

```bash
java -jar target/transept-core-0.1-SNAPSHOT-jar-with-dependencies.jar repair \
  --topology ../Dataset_Structural_GPCRs/GPCRdb/sequences/tmregions/classF_predicted_topologies.3line \
  --output-dir repaired \
  ../Dataset_Structural_GPCRs/GPCRdb/precomputed/classF/mafft_classF.fasta
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--topology` | (required) | Three-line topology file of the proteins of the alignments |
| `--output-dir` | (required) | Directory for the repaired alignments; it must not be that of an input |
| `ALIGNMENT...` | (required) | One or more alignments (FASTA) of the proteins of the topology |

An alignment whose rows are not the proteins of the topology, with their residues, is an error.
