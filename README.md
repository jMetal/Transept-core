# TM-M2StructAlign: Multiobjective Structure-Guided Alignment of Transmembrane Proteins

[![Dataset DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23243974.svg)](https://doi.org/10.5281/zenodo.23243974)

_A multiobjective evolutionary algorithm for multiple sequence alignment (MSA) of G
protein-coupled receptors guided by predicted topology and AlphaFold2 structures._

## Overview

TM-M2StructAlign refines alignments produced by existing tools (MAFFT, Kalign, ClustalW...)
with a multiobjective evolutionary algorithm built on [jMetal](https://github.com/jMetal/jMetal):
NSGA-II in the published configuration, and other algorithms for comparison. Every solution is
an alignment, and two objectives are maximised at the same time:

| Objective | What it rewards |
| --- | --- |
| **SOPwTP** (topology-aware sum of pairs) | Residue pairs scored with PHAT when both residues are transmembrane (TM) and BLOSUM62 otherwise; affine gap penalties that are higher in TM regions (8/3) than elsewhere (3/1) |
| **lDDT** (structural consistency) | Superposition-free lDDT computed from AlphaFold2 Cα distance matrices (inclusion radius 15 Å): aligned residues should have the same local distances in every structure |

The search never cuts a TM helix when recombining alignments, and the default mutation moves
gap blocks preferably into loops. The result is a set of Pareto-optimal alignments that trade sequence and structural
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

## Getting started

Requirements: Java 21 and Maven 3.9+. Dependencies (jMetal, picocli, BioJava) are downloaded from
Maven Central; until jMetal 7.11 is released, install `7.11-SNAPSHOT` locally from the jMetal
sources (`mvn install -DskipTests`).

```bash
git clone https://github.com/TM-M2StructAlign/TM-M2StructAlign.git
git clone https://github.com/TM-M2StructAlign/Dataset_Structural_GPCRs.git
git -C Dataset_Structural_GPCRs checkout v1.0.0
cd TM-M2StructAlign
mvn clean package -DskipTests      # target/TM-MSAStructAligner-1.0-jar-with-dependencies.jar

java -jar target/TM-MSAStructAligner-1.0-jar-with-dependencies.jar run \
  --data-root ../Dataset_Structural_GPCRs/GPCRdb \
  --dataset classF \
  --seeds clustalw,kalign,mafft,muscle5,tcoffee,tmaligner \
  --output-root results --run-id first --seed 1 \
  --max-evaluations 500 --population-size 20
```

The run writes the alignments of the final population (`MSASol<i>.fasta`) and their objective
values (`FUN.csv`) to `results/ejecuciones/classF/first/`. [T1](docs/tutorials/T1-first-run.md)
explains every step and file.

## Documentation

| Tutorial | What you will learn |
| --- | --- |
| [T1. First run in five minutes](docs/tutorials/T1-first-run.md) | Build, run on the smallest dataset, read the results, repeat a run exactly |
| [T2. Problem formulation, encoding and variation operators](docs/tutorials/T2-problem-formulation.md) | The objectives, the constraints, how alignments are encoded, and what every crossover and mutation does, with a runnable example |
| T3. Reproducing the published experiments | *In preparation* |
| [T4. Choosing an alignment and evaluating it against GPCRdb](docs/tutorials/T4-choosing-and-evaluating.md) | Reading the non-dominated alignments, extremes and knee, and comparing them and the seeds with the GPCRdb reference alignments |
| T5. Comparing algorithms and operators | *In preparation* |
| [T6. Using your own dataset](docs/tutorials/T6-your-own-dataset.md) | The input files a dataset needs, how to compute the distance matrices, and how to check that structures and seeds are read correctly |
| T7. Adding a variation operator | *In preparation* |

Reference:

- [Command line](docs/reference/cli.md): every option of `run` and `baseline`, the algorithms and
  operators, and the comparison with reference alignments.
- [Output files](docs/reference/outputs.md): the files of a run and their formats.
- [Metrics](docs/reference/metrics.md): the objectives and the reference-based metrics.

## Development

- Read [AGENTS.md](AGENTS.md), [JAVA_CODING_GUIDELINES.md](JAVA_CODING_GUIDELINES.md) and
  [GIT_GUIDELINES.md](GIT_GUIDELINES.md).
- `mvn verify` enforces Google Java Style (Checkstyle) and 80 % line coverage on
  `org.tmm2structalign`.
- Integration tests that need the real dataset run only when `GPCR_DATASET` points to the
  `GPCRdb` directory. `BaselineOracleIT` checks that the objectives reproduce the values of the
  original implementation stored in the dataset's `baselines.csv`.
- User documentation lives in `docs/tutorials/` and `docs/reference/`, in Markdown.
  `FormulationTourTest` fails if tutorial T2 no longer shows the output of
  `example.FormulationTour`: update the tutorial when an objective or an operator changes.
- The state before the November 2026 rewrite (Java 17, jMetal 6.2, original runner) is tagged
  `legacy-v1.0`.

## License

Licensed under the terms in [LICENSE](LICENSE).
