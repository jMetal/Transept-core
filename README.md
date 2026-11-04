# TM-M2StructAlign: Multiobjective Structure-Guided Alignment of Transmembrane Proteins

[![Dataset DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23243974.svg)](https://doi.org/10.5281/zenodo.23243974)

_A multiobjective evolutionary algorithm for multiple sequence alignment (MSA) of G
protein-coupled receptors guided by predicted topology and AlphaFold2 structures._

## Overview

TM-M2StructAlign refines alignments produced by existing tools (MAFFT, Kalign, ClustalW...)
with an NSGA-II built on [jMetal](https://github.com/jMetal/jMetal). Every solution is an
alignment, and two objectives are optimised at the same time:

| Objective | What it rewards |
| --- | --- |
| **SOPwTP** (topology-aware sum of pairs) | Residue pairs scored with PHAT when both residues are transmembrane (TM) and BLOSUM62 otherwise; affine gap penalties that are higher in TM regions (8/3) than elsewhere (3/1) |
| **lDDT** (structural consistency) | Superposition-free lDDT computed from AlphaFold2 Cα distance matrices (inclusion radius 15 Å): aligned residues should have the same local distances in every structure |

The search never cuts a TM helix when recombining alignments, and moves gap blocks preferably
into loops. The result is a set of Pareto-optimal alignments that trade sequence and structural
consistency.

## Project ecosystem

TM-M2StructAlign is maintained together with the benchmark dataset and the manuscript sources:

| Resource | Repository |
| --- | --- |
| **TM-M2StructAlign source code** | [TM-M2StructAlign](https://github.com/TM-M2StructAlign/TM-M2StructAlign) |
| **GPCR structural benchmark and experimental resources** | [Dataset_Structural_GPCRs](https://github.com/TM-M2StructAlign/Dataset_Structural_GPCRs) · [Zenodo v1.0.0](https://doi.org/10.5281/zenodo.23243974) |
| **Main manuscript and supplementary material** | [TM-M2StructAlign-Paper](https://github.com/JOELITO07/TM-M2StructAlign-Paper) |

### Related publications

1. Cedeño-Muñoz, J., Zambrano-Vega, C., and Nebro, A. J. (2025). **TMP-M2Align: A Topology-Aware Multiobjective Approach to the Multiple Sequence Alignment of Transmembrane Proteins.** *Algorithms*, 18(10), 640. [https://doi.org/10.3390/a18100640](https://doi.org/10.3390/a18100640)
2. Cedeño-Muñoz, J., Zambrano-Vega, C., and Nebro, A. J. **TM-M2StructAlign: A Multiobjective Tool for Structure-Guided Multiple Sequence Alignment of G Protein-Coupled Receptors Using AlphaFold2-Derived Constraints.** Manuscript prepared for *Computational Biology and Chemistry*. Source and supplementary material: [TM-M2StructAlign-Paper](https://github.com/JOELITO07/TM-M2StructAlign-Paper).
3. Zambrano-Vega, C., Cedeño-Muñoz, J., and Nebro, A. J. (2026). **A Curated Dataset of Human GPCR Sequences, AlphaFold-Predicted Structures, Transmembrane Topologies, and Reference Alignments for Structural Bioinformatics** (Version 1.0.0) [Dataset]. Zenodo. [https://doi.org/10.5281/zenodo.23243974](https://doi.org/10.5281/zenodo.23243974).

A companion data article documenting the benchmark resource is maintained in [TM-M2StructAlign_Dataset_paper](https://github.com/JOELITO07/TM-M2StructAlign_Dataset_paper).

### Versioned benchmark dataset

The benchmark dataset associated with TM-M2StructAlign is versioned as **v1.0.0** and assigned the version-specific Zenodo DOI **[10.5281/zenodo.23243974](https://doi.org/10.5281/zenodo.23243974)**. The corresponding GitHub release is [Dataset_Structural_GPCRs v1.0.0](https://github.com/TM-M2StructAlign/Dataset_Structural_GPCRs/releases/tag/v1.0.0), based on source snapshot `75287f1da907248aa29b8fb5adf17a4c63891797`.

For reproducible analyses, cite the Zenodo version-specific DOI and use the archived v1.0.0 dataset snapshot rather than a later state of the live GitHub repository.

---

## How it works

- **Representation.** One gap layout per protein: the sorted runs of gap columns of its aligned
  row (`alignment.GapIntervals`). Residues are never stored; the alignment is decoded from the
  ungapped sequences.
- **Initial population.** An explicit list of seed alignments (precomputed with other tools),
  without duplicates, completed with crossover children up to the population size.
- **Crossover** (`TmSafeSinglePointCrossover`). Single cut at a column of the first parent; the
  cut is rejected if it separates two consecutive TM residues of any protein. Remaining residues
  of every protein are taken from the second parent.
- **Mutation** (`TmAwareGapShiftMutation`). Moves a gap block of one protein to a new position
  chosen by rejection sampling: TM residues are avoided (acceptance 0.05), loops and residues
  without structure are preferred (0.90 / 0.95).
- **Structure mapping.** AlphaFold residues are mapped to the input sequence by PDB numbering
  when the residues agree, or by sequence alignment otherwise (e.g. OXER1 and MCHR1, whose
  AlphaFold model is a different isoform). Every mapping is reported in the run metadata.

## Requirements and build

- Java 21 and Maven 3.9+.
- Dependencies (jMetal 7.11, picocli, BioJava) are downloaded from Maven Central. Until jMetal
  7.11 is released, install `7.11-SNAPSHOT` locally from the jMetal sources (`mvn install`).

```bash
mvn clean verify     # compile, unit and integration tests, Checkstyle, coverage
mvn clean package    # target/TM-MSAStructAligner-1.0-jar-with-dependencies.jar
```

The benchmark data lives in the companion repository
[Dataset_Structural_GPCRs](https://github.com/TM-M2StructAlign/Dataset_Structural_GPCRs); the
commands below take its `GPCRdb` directory as `--data-root`.

## Usage

### Run the algorithm

```bash
java -jar target/TM-MSAStructAligner-1.0-jar-with-dependencies.jar run \
  --data-root ../Dataset_Structural_GPCRs/GPCRdb \
  --dataset classF \
  --seeds clustalw,kalign,mafft,muscle5,tcoffee,tmaligner \
  --output-root results \
  --run-id seed_000 \
  --seed 0 \
  --max-evaluations 10000 \
  --population-size 100 \
  --threads 4
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--seeds` | (required) | Tools whose `precomputed/<dataset>/<tool>_<dataset>.fasta` seed the population |
| `--seed` | (required) | Random seed; the same seed and options reproduce the run exactly |
| `--max-evaluations` | 10000 | Stopping criterion |
| `--population-size` | 100 | Population and offspring size |
| `--threads` | 1 | Parallel evaluation threads (results do not depend on it) |
| `--algorithm` | `NSGA_II` | Search algorithm (see below) |
| `--initialization` | `CROSSOVER` | Non-seed members: `CROSSOVER` (original), `CROSSOVER_DISTINCT` (no duplicates after removing gap-only columns) or `PERTURBED` (crossover child + 10 mutations) |
| `--crossover` | `SINGLE_POINT` | `SINGLE_POINT` (original) or `TWO_POINT`; both never cut a TM helix |
| `--mutation` | `GAP_SHIFT` | `GAP_SHIFT` (original: discourages gaps inside TM helices). The rest never break a helix: `TM_SAFE_GAP_SHIFT`, `INSERT_GAP` (one gap per row, +1 column), `SPLIT_RESIDUE_RUN` (new one-gap run between residues of every row, +1 column), `MERGE_GAP_RUNS` (joins two consecutive runs of a row) and `RANDOM_CHOICE` (one of the four, uniformly) |
| `--crossover-probability` / `--mutation-probability` | 0.8 / 0.2 (1.0 for PAES) | Operator probabilities |
| `--mapping-policy` | `VERIFIED` | `VERIFIED`, `NUMBERED_ONLY` or `LEGACY` (original, unchecked mapping) |
| `--lddt-aggregation` | `HARMONIC_OVER_LENGTH` | Original aggregation, or `PAIR_MEAN` |
| `--progress-interval` | 100 | Evaluations between progress records |

Options can be stored in a file, one per line, and passed as `@run-classF.args`. A batch of ten
seeds on Windows is available in `scripts/run_10_seeds_windows.ps1`.

### Algorithms

All algorithms use the same objectives, operators and seeded initial population; they stop after
`--max-evaluations`.

| `--algorithm` | Algorithm | Notes |
| --- | --- | --- |
| `NSGA_II` | NSGA-II | Published configuration (default) |
| `NSGA_III` | NSGA-III | Population size must be a multiple of 4 |
| `AGE_MOEA`, `AGE_MOEA_II` | AGE-MOEA, AGE-MOEA-II | |
| `SMS_EMOA` | SMS-EMOA | Steady state: `--threads` does not speed it up |
| `MOEA_D` | MOEA/D | Normalised PBI aggregation; steady state |
| `PAES` | PAES, (1+1) with a crowding-distance archive | Mutation only; starts from the first seed; archive size = population size |
| `RANDOM_SEARCH` | Random search, gaps at random | Sanity check (weak baseline) |
| `RANDOM_SEEDED` | Random search, seeds + 10 random mutations | Sanity check (strong baseline) |

SMS-EMOA, MOEA/D and PAES need jMetal 7.11, which lets them start from a given population or
solution.

### Score precomputed alignments

```bash
java -jar target/TM-MSAStructAligner-1.0-jar-with-dependencies.jar baseline \
  --data-root ../Dataset_Structural_GPCRs/GPCRdb \
  --datasets classF,classT2_19 \
  --output baselines.csv
```

### Compare with the GPCRdb reference alignments

`org.tmm2structalign.evaluation.EvaluationMain` computes Pair F1, TM-Pair F1 and TM-gap rate
for every alignment of a result tree; see [docs/evaluation-metrics.md](docs/evaluation-metrics.md).
`scripts/pareto_analysis.py` analyses the fronts of several runs.

## Output files

Each run writes to `<output-root>/ejecuciones/<dataset>/<run-id>/` and never overwrites an
existing directory.

| File | Content |
| --- | --- |
| `FUN.csv` | One line per solution: SOPwTP and lDDT, comma-separated, no header (jMetal convention) |
| `MSASol<i>.fasta` | Alignment of line `i` of `FUN.csv` |
| `config.json` | Full run configuration |
| `metadata.json` | Git commit, timing, seed files with SHA-256 (and duplicates), structure mapping of every protein, proteins without structure |
| `progress.jsonl` | One JSON line per progress record (evaluations, elapsed time, best values) |
| `runtime.csv` | `runtime_ms,runtime_seconds,evaluations` with one row |
| `seed.txt` | Random seed |
| `error.json` | Only if the run failed: exception type and message |

## Project layout

```text
src/main/java/org/tmm2structalign/
  model/       proteins, topology
  alignment/   gap layouts, decoded alignments, jMetal solution
  io/          3-line topology, FASTA, distance matrices, PDB readers
  structure/   sequence-structure mapping and structural neighbourhoods
  objective/   SOPwTP, lDDT, substitution matrices
  operator/    crossover and mutation
  problem/     jMetal problem
  algorithm/   seeded population and NSGA-II assembly
  run/         run configuration, execution, outputs, baseline scoring
  cli/         command line (picocli)
  evaluation/  comparison with reference alignments
src/main/java/org/tm_msaligner/   legacy irace tuning code (kept as is)
docs/review/                       code and dataset review, work plan
```

## Development

- Read [AGENTS.md](AGENTS.md), [JAVA_CODING_GUIDELINES.md](JAVA_CODING_GUIDELINES.md) and
  [GIT_GUIDELINES.md](GIT_GUIDELINES.md).
- `mvn verify` enforces Google Java Style (Checkstyle) and 80 % line coverage on
  `org.tmm2structalign`.
- Integration tests that need the real dataset run only when `GPCR_DATASET` points to the
  `GPCRdb` directory. `BaselineOracleIT` checks that the objectives reproduce the values of the
  original implementation stored in the dataset's `baselines.csv`.
- The state before the November 2026 rewrite (Java 17, jMetal 6.2, original runner) is tagged
  `legacy-v1.0`.

## License

Licensed under the terms in [LICENSE](LICENSE).
