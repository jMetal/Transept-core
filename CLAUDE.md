# CLAUDE.md

Transept core is the engine of [Transept](https://github.com/jMetal/Transept): a Java 21 command
line (picocli, jMetal 7.11, Evolver 2.6) that aligns transmembrane proteins with multiobjective
evolutionary algorithms, guided by their topology and AlphaFold2 structures. The Transept graphical
interface (Streamlit, separate repository) drives it only through its commands and output files,
so those are a contract between the two projects.

Both projects are developed collaboratively and are still prototypes (0.1): options, output
formats and defaults may change, but every change must keep the two repositories in step.

## Current state

- **Commands**: `run`, `seeds`, `algorithms`, `baseline`, `score` and `repair`
  (`docs/reference/cli.md`).
- **Formulations** (`--objectives`): SOPwTP and lDDT (default); SOPwTP, lDDT and aligned-pair
  coverage (experimental); lDDT alone for the genetic algorithms.
- **Algorithms**: NSGA-II (default), SPEA2 and a GA built with Evolver from parameter spaces; their
  asynchronous variants; NSGA-III, AGE-MOEA, AGE-MOEA-II, SMS-EMOA, MOEA/D, PAES and random searches.
- **Operators** never split a transmembrane helix. Gaps inside helices that come from the seed
  alignments are moved out first by default; `run --no-repair-seeds` keeps those of the precomputed
  alignments.
- **Seeds**: precomputed alignments of the dataset (`--seeds`) or internal ones computed with BioJava.
- **Evaluation**: `score` compares alignments with a reference (Pair F1, TM-pair F1, TM-gap rate,
  coverage). The benchmark is Dataset_Structural_GPCRs v1.0.0 (Zenodo), human GPCRs with GPCRdb
  reference alignments.
- **Quality**: `mvn clean verify` runs unit and integration tests, Checkstyle and JaCoCo (80 %);
  GitHub Actions runs it on Linux.
- **Releases**: none yet. The tag `legacy-v1.0` keeps the state of TM-M2StructAlign, from which this
  project derives. See `CHANGELOG.md` (`[Unreleased]`).

## Open lines of work

Contributions are welcome in any of these; open an issue first to agree on the approach.

- A first release: versioned jar, archived with a DOI, and an installation that works on a clean
  macOS or Linux machine.
- Transmembrane families beyond GPCRs (e.g. the transmembrane reference set of BAliBASE), with
  topology and AlphaFold models from UniProt.
- Seeds from structural aligners (FoldMason, US-align, MUSTANG) and their effect on the results.
- Known issues of the benchmark data: the GPCRdb reference of classB2 truncates `agrv1_human`; the
  structural tools' outputs of classA_002/003/004 do not match the GPCRdb sequences. Fixes go to a
  new version of the dataset, never to v1.0.0.
- New operators or algorithms, evaluated against the defaults (tutorials T5 and T7).

## Getting started

```bash
# Side by side: this repository, the interface (optional) and the benchmark data
git clone https://github.com/jMetal/Transept-core.git
git clone https://github.com/jMetal/Transept.git
git clone https://github.com/TM-M2StructAlign/Dataset_Structural_GPCRs.git
git -C Dataset_Structural_GPCRs checkout v1.0.0

cd Transept-core
git switch develop
mvn clean verify                                              # must pass before any commit
GPCR_DATASET=../Dataset_Structural_GPCRs/GPCRdb mvn verify    # also the tests on real data
```

A run of a few seconds: `classF` with 500 evaluations and a population of 20 (skill `quick-run`,
or tutorial `docs/tutorials/T1-first-run.md`).

## How to contribute

- Start a branch from `develop` (`feat/...`, `fix/...`, `docs/...`) and open a pull request into
  `develop`. Only the maintainers merge `develop` into `main`.
- Commits follow `GIT_GUIDELINES.md` (Conventional Commits, one logical change per commit, code,
  tests and docs in separate commits, every commit builds).
- A change to a command, an option, an output file or a parameter space is a change of the contract
  with the interface: update the documentation here and the interface in a pull request of its own
  (skill `cli-contract`).
- Changes to the objectives, the operators' behaviour or the default configuration change the
  results: propose them in an issue or in the pull request, with the evidence, before making them.

Project skills (`.claude/skills/`): `check` runs every check before a commit or a pull request;
`cli-contract` lists what to update when the command line or the outputs change; `quick-run`
builds the jar and runs a short example end to end.

## Rules for the assistant

- Never push, open pull requests, create tags or publish anything unless the user asks, each time.
  Never commit to `main`.
- Run `mvn clean verify` before proposing a commit, and report failures with their output.
- Do not commit large data or run outputs, and never modify the published v1.0.0 files of
  Dataset_Structural_GPCRs.
- Code, comments, commit messages and documentation in English.

The rules, the code map and the pitfalls of the project are in the files below; read them before
changing code.

@AGENTS.md
@GIT_GUIDELINES.md
@JAVA_CODING_GUIDELINES.md
