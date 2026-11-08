# Transept core

_Multiobjective, structure-guided multiple sequence alignment of transmembrane proteins._

Transept core is the engine of [Transept](https://github.com/jMetal/Transept), a tool for
aligning transmembrane proteins such as G protein-coupled receptors (GPCRs). It is a Java command
line that can be used on its own or driven by the Transept graphical interface.

> **Status: prototype (0.1.0).** Options, output formats and defaults may still change. Please
> report problems and suggestions in the [issues](https://github.com/jMetal/Transept-core/issues).

## What it does

Transept refines alignments produced by existing tools (MAFFT, Kalign, ClustalW...) or computed
internally with BioJava, using multiobjective evolutionary algorithms built on
[jMetal](https://github.com/jMetal/jMetal) and [Evolver](https://github.com/jMetal/Evolver).
Every solution is an alignment, and by default two objectives are maximised at the same time:

| Objective | What it rewards |
| --- | --- |
| **SOPwTP** (topology-aware sum of pairs) | Residue pairs scored with PHAT when both residues are transmembrane (TM) and BLOSUM62 otherwise; affine gap penalties that are higher in TM regions (8/3) than elsewhere (3/1) |
| **lDDT** (structural consistency) | Superposition-free lDDT computed from AlphaFold2 Cα distance matrices (inclusion radius 15 Å): aligned residues should have the same local distances in every structure |

`--objectives` selects other formulations (the lDDT alone, or a third objective, the
aligned-pair coverage). The search never cuts a TM helix when recombining alignments, and the
default mutation moves gap blocks preferably into loops. The result is a set of Pareto-optimal
alignments that trade sequence and structural consistency.

Algorithms: NSGA-II (default), SPEA2 and a genetic algorithm, configurable through Evolver
parameter spaces, their asynchronous parallel variants, and NSGA-III, AGE-MOEA, SMS-EMOA,
MOEA/D, PAES and random search for comparison (`java -jar ... algorithms` lists them).

## Getting started

Requirements: Java 21 or later.

The simplest way is to download `transept-core-0.1.0-jar-with-dependencies.jar` from the
[releases](https://github.com/jMetal/Transept-core/releases) and the benchmark data:

```bash
git clone https://github.com/TM-M2StructAlign/Dataset_Structural_GPCRs.git
git -C Dataset_Structural_GPCRs checkout v1.0.0

java -jar transept-core-0.1.0-jar-with-dependencies.jar run \
  --data-root Dataset_Structural_GPCRs/GPCRdb \
  --dataset classF \
  --seeds clustalw,kalign,mafft,muscle5,tcoffee,tmaligner \
  --output-root results --run-id first --seed 1 \
  --max-evaluations 500 --population-size 20
```

The run writes the alignments of the final population (`MSASol<i>.fasta`) and their objective
values (`FUN.csv`) to `results/ejecuciones/classF/first/`. [T1](docs/tutorials/T1-first-run.md)
explains every step and file.

### Building from source

Requirements: Java 21 and Maven 3.9+. Every dependency (jMetal 7.11, Evolver 2.6, BioJava,
picocli) is downloaded from Maven Central.

```bash
git clone https://github.com/jMetal/Transept-core.git && cd Transept-core
mvn clean package -DskipTests      # target/transept-core-0.1.0-jar-with-dependencies.jar
```

## Benchmark data

The examples and tutorials use the GPCR structural benchmark
[Dataset_Structural_GPCRs](https://github.com/TM-M2StructAlign/Dataset_Structural_GPCRs),
version **v1.0.0** ([Zenodo DOI 10.5281/zenodo.23243974](https://doi.org/10.5281/zenodo.23243974)).
For reproducible analyses, use that archived version rather than a later state of the live
repository. [T6](docs/tutorials/T6-your-own-dataset.md) explains how to prepare your own data.

## Documentation

| Tutorial | What you will learn |
| --- | --- |
| [T1. First run in five minutes](docs/tutorials/T1-first-run.md) | Build, run on the smallest dataset, read the results, repeat a run exactly |
| [T2. Problem formulation, encoding and variation operators](docs/tutorials/T2-problem-formulation.md) | The objectives, the constraints, how alignments are encoded, and what every crossover and mutation does |
| [T3. A study with several seeds](docs/tutorials/T3-several-seeds.md) | Running several seeds, combining their fronts and choosing alignments from the combined front |
| [T4. Choosing an alignment and evaluating it against GPCRdb](docs/tutorials/T4-choosing-and-evaluating.md) | Reading the non-dominated alignments and comparing them with the GPCRdb reference alignments |
| [T5. Comparing algorithms and operators](docs/tutorials/T5-comparing-algorithms-and-operators.md) | Running configurations fairly, measuring them with quality indicators and paired tests |
| [T6. Using your own dataset](docs/tutorials/T6-your-own-dataset.md) | The input files a dataset needs and how to check that they are read correctly |
| [T7. Adding a variation operator](docs/tutorials/T7-adding-an-operator.md) | The rules every operator follows, its tests, and how to expose it in the command line |

Reference: [command line](docs/reference/cli.md),
[configuring the algorithms](docs/reference/configuration.md),
[operators](docs/reference/operators.md), [output files](docs/reference/outputs.md) and
[metrics](docs/reference/metrics.md).

## Related publications

1. Cedeño-Muñoz, J., Zambrano-Vega, C., and Nebro, A. J. (2025). **TMP-M2Align: A Topology-Aware
   Multiobjective Approach to the Multiple Sequence Alignment of Transmembrane Proteins.**
   _Algorithms_, 18(10), 640. <https://doi.org/10.3390/a18100640>
2. Zambrano-Vega, C., Cedeño-Muñoz, J., and Nebro, A. J. (2026). **A Curated Dataset of Human GPCR
   Sequences, AlphaFold-Predicted Structures, Transmembrane Topologies, and Reference Alignments
   for Structural Bioinformatics** (Version 1.0.0) [Dataset]. Zenodo.
   <https://doi.org/10.5281/zenodo.23243974>

## Development

- Read [AGENTS.md](AGENTS.md), [JAVA_CODING_GUIDELINES.md](JAVA_CODING_GUIDELINES.md) and
  [GIT_GUIDELINES.md](GIT_GUIDELINES.md).
- `mvn verify` enforces Google Java Style (Checkstyle) and 80 % line coverage.
- Integration tests that need the real dataset run only when `GPCR_DATASET` points to the
  `GPCRdb` directory.
- `scripts/slurm/` prepares and runs experimental campaigns on a Slurm cluster.

## License

Transept core is free software, distributed under the
[GNU General Public License v3.0](LICENSE).
