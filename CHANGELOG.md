# Changelog

All notable changes to Transept core are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.9] - 2026-11-27

First release.

Transept core succeeds TM-M2StructAlign. The state of that tool before the changes below is kept
in the tag `legacy-v1.0`, which reproduces its earlier results.

### Added

- Publication in Maven Central as `org.uma.jmetal:transept-core`. The jar with dependencies is
  built locally and not deployed.
- Command line with the commands `run`, `seeds`, `algorithms`, `baseline`, `score` and `repair`.
- `repair`: moves the gaps that split a transmembrane helix out of any alignment, without search,
  as `run` does with the seeds (tutorial T8).
- `run` repairs the seeds by default; `--no-repair-seeds` keeps the gaps of the precomputed ones. The
  internal seeds are always repaired.
- Algorithms: NSGA-II (default), SPEA2 and a genetic algorithm configurable from Evolver parameter
  spaces; their asynchronous parallel variants; NSGA-III, AGE-MOEA, AGE-MOEA-II, SMS-EMOA, MOEA/D,
  PAES, and two random searches as sanity checks.
- The offspring size of the configurable algorithms is 100 by default, whatever the population size.
  NSGA-III and MOEA/D need weight vectors for their population size with three objectives.
- Formulations: SOPwTP and lDDT (default), lDDT alone for the genetic algorithms, and SOPwTP, lDDT
  and aligned-pair coverage (experimental).
- Operators that never split a transmembrane helix: single-point and uniform crossovers; gap shift,
  merge of gap runs, block shift and window realignment mutations, and a random choice among these
  four, with their parameters in the configurations.
- Internal seed alignments computed with BioJava, so that a dataset needs no precomputed alignments.
- `baseline` scores seed alignments, including those given with `--alignments`; `score` compares
  alignments with a reference alignment (pair F1, TM-pair F1, TM-gap rate, coverage).
- Records of every run: `config.json`, `metadata.json` (build, seeds with SHA-256, mapping status of
  every structure), `progress.jsonl` and `error.json`.
- A run can be stopped early, keeping its results, by creating a `STOP` file in its directory.
- Documentation: README, seven tutorials and references of the command line, configurations,
  operators, metrics and outputs.

### Changed

- The structural objective is the mean of the per-residue lDDT scores; the original aggregation
  favoured short, compact alignments.
- Seed alignments are given explicitly with `--seeds` and recorded, instead of reading every file of
  a folder.
- The initial population holds distinct seeds and distinct crossover children, of the exact size.
- Java 21 and jMetal 7.11; Google Java Style and 80 % line coverage are checked by `mvn verify`.

### Fixed

- Invalid input of `repair`, `baseline`, `score` and `seeds` (an alignment of other proteins, a
  missing file) printed a Java stack trace; it is now one line, `Error: <message>`, with exit code 1.
- The AlphaFold distances are mapped to the sequence checking residue identity, not only the
  numbering: isoforms that differ from the sequence were mapped to wrong residues.
- The window realignment scores residue pairs with PHAT only when both are transmembrane, as SOPwTP
  does; one was enough, and an ambiguous residue (X, B, Z) facing a helix made the run fail.
- Gap-only columns are always removed before evaluation.
- The gap-shift mutation could loop forever when no row had a run of two or more gaps.
- Proteins without a structure are recorded in `metadata.json` instead of only being reported on
  the console.
