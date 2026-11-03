# AGENTS.md

Guidance for AI coding assistants (and humans) working on **TM-M2StructAlign**. Read this file
first, then `JAVA_CODING_GUIDELINES.md` and `GIT_GUIDELINES.md`.

## 1. What this project is

TM-M2StructAlign is a multiobjective evolutionary algorithm (NSGA-II-style, built on jMetal) for
multiple sequence alignment (MSA) of transmembrane proteins, specifically human GPCRs. It
optimizes two objectives:

1. **SOPwTP** (`score/impl/StructSumOfPairsWithTopologyPredict`): sum of pairs using PHAT inside
   TM regions and BLOSUM62 elsewhere, with affine gap penalties that differ for TM and non-TM.
2. **Structural lDDT** (`score/impl/LDDTStructuralScore`): reference-free, superposition-free
   lDDT computed from AlphaFold2 Cα distance matrices (R0 = 15 Å).

Both are maximized; internally they are negated for jMetal minimization.

## 2. Workspace layout

Two sibling repositories are expected under the same parent directory:

```text
MSA/
├── TM-M2StructAlign/            # this repository (Java, Maven)
└── Dataset_Structural_GPCRs/    # benchmark data (separate git repo, Zenodo v1.0.0)
    └── GPCRdb/
        ├── sequences/{fasta,tmregions,gff3}/   # <ds>.fasta, <ds>_predicted_topologies.3line
        ├── distances/<ds>/<prot>_D.csv, <prot>_idx.csv
        ├── precomputed/<ds>/*.fasta            # seed alignments + baselines.csv
        ├── reference_alignments/<ds>.{fasta,msf}
        ├── alphafold_pdb/, alphafold_json/, weights/
        └── ...
```

Datasets (`<ds>`): `classA_001..004`, `classB1`, `classB2`, `classC`, `classF`, `classT2`, and
the reduced `_19` variants (19 proteins each, required by mTM-align).

## 3. Build, test, run

```bash
mvn clean test                     # must pass before every commit
mvn clean package                  # builds target/TM-MSAStructAligner-1.0-jar-with-dependencies.jar
```

Main runner (`org.tm_msaligner.TM_MSAStructAlignerGPCRdb`, 9 positional arguments):

```bash
java -jar target/TM-MSAStructAligner-1.0-jar-with-dependencies.jar \
  <dataDirectory=.../GPCRdb> <dataset> <maxEvaluations> <populationSize> \
  <numberOfCores> <observerFrequency> <runId> <outputRoot> <seed>
```

Output goes to `<outputRoot>/ejecuciones/<dataset>/<runId>/` (`FUN.tsv`, `MSASol*.fasta`,
`BestScores.tsv`, `seed.txt`, `runtime.txt`). The runner refuses to overwrite an existing folder.

Other entry points: `BaselineMain` (scores every precomputed FASTA; currently uses a hard-coded
Windows path), `evaluation/EvaluationMain` (Pair F1, TM-Pair F1, TM-gap rate against GPCRdb
references; see `docs/evaluation-metrics.md`), `scripts/pareto_analysis.py`.

For quick experiments use a small dataset (`classF`, 11 sequences) with ~400 evaluations and a
population of 20: it runs in seconds.

## 4. Code map

| Package | Role |
| --- | --- |
| `problem`, `problem.impl` | Problem loading (3-line topology, distance/idx matrices, seed alignments) and evaluation. Structural variant: `StructuralTMMSAProblem`, `MultiObjStructuralTMMSAProblem` |
| `solution` | Gap-group encoding of MSAs (`StructuralTM_MSASolution`): per sequence, a sorted list of `[start, end]` gap intervals |
| `crossover` | `BioSPXMSACrossover`: single-point crossover that never cuts inside a TM helix |
| `mutation` | `BioShiftClosedGapsMSAMutation`: moves a gap block, biased away from TM residues |
| `score.impl` | Objective functions |
| `algorithm.structural_multiobjective` | `StructuralTM_M2Align` and its builder (jMetal component-based EA) |
| `solutionscreation` | Initial population from precomputed alignments + crossover |
| `evaluation` | Reference-based metrics |
| `auto`, `parallel`, `algorithm.{single,multi}objective` | Earlier, non-structural TM-M2Align variants and irace tuning |
| `test` (under `src/main`!) | Legacy experiment drivers with hard-coded paths; do not extend |

## 5. Conventions

- Follow `JAVA_CODING_GUIDELINES.md` for **new and modified code**. Platform: Java 21, jMetal
  7.10, JUnit 6 (Jupiter) + Mockito (`mock`/`when`/`verify` only), unit tests `*Test` (surefire)
  and integration tests `*IT` (failsafe, `mvn verify`).
- The rewritten core lives in **`org.tmm2structalign`**. Checkstyle (Google style) and JaCoCo
  (>= 80 % line coverage) are enforced on that package only; `mvn verify` fails otherwise.
  The legacy package `org.tm_msaligner` is exempt and is being replaced (see §7).
- jMetal 7.10 is not yet on Maven Central (2026-11-03); builds resolve it from `~/.m2`. Check that
  it resolves from Central before pushing.
- Follow `GIT_GUIDELINES.md`: Conventional Commits, one logical change per commit, never mix
  `test:` with `fix:`/`feat:`, never mix code with docs.
- Write code, comments, and messages in English. Discussion and review documents may be in
  Spanish.
- Use `JMetalRandom` for every stochastic decision so that runs are reproducible from the seed.

## 6. Pitfalls (read before touching these areas)

1. **Sequence→structure mapping.** `<prot>_idx.csv` holds AlphaFold PDB residue numbers, which
   the loader treats as 1-based positions in the GPCRdb sequence **without checking residue
   identity**. That is wrong for OXER1 and MCHR1, whose AlphaFold sequence differs. See A1.
2. **`removeGapColumns` is buggy**: it leaves gap-only columns (reproduced for 2 sequences and for
   ≥3 sequences with several gap blocks). Because the lDDT divides by alignment length, this
   changes objective 2. See A2.
3. **Seed alignments**: every `*.fasta` in `precomputed/<ds>/` is loaded, and one invalid file
   aborts the run. Currently classA_002/003/004 abort because of the structural-tool outputs.
   Caretta outputs are residue-level misaligned for most classes. See A4 and
   `docs/review/01-caretta-mtmalign.md`.
4. **The lDDT objective rewards compact alignments** and is nearly blind to helix register shifts.
   Do not draw conclusions from objective 2 alone; check against reference-based metrics.
5. `.gitignore` excludes `*.fasta` except under `src/test/resources/`. Prefer generating fixtures
   in `@TempDir`; version only small, hand-checked files.
6. Missing structural files only produce a `[WARN]` on stdout; the sequence is silently excluded
   from the lDDT.

## 7. Current work (November 2026)

The goal on branch `antonio` is to (1) review TM-M2StructAlign for errors and inconsistencies,
(2) review the dataset directory, and (3) carry out an improvement plan. Documents:

- `docs/review/01-caretta-mtmalign.md`: investigation of anomalous baselines
- `docs/review/02-code-review.md`: code findings (ids A*, M*, B*), plus the JDK 21 / jMetal 7.9
  check
- `docs/review/03-dataset-review.md`: dataset findings (ids D*)
- `docs/review/04-work-plan.md`: phased plan; **Phase 0 decisions belong to the authors**. Do
  not implement work that depends on an open decision without asking.
- `docs/review/scripts/`: read-only diagnostic scripts used to produce the numbers

An approved rewrite is under way (November 2026): the structural core is re-implemented with
tests in `org.tmm2structalign` (tag `legacy-v1.0` marks the pre-rewrite state), then the old
structural flow in `org.tm_msaligner` is removed. The irace code and its dependencies stay (§8).

When fixing a finding, reference its id in the commit body (e.g. `Refs: A2`) and update the
corresponding document.

## 8. Do not

- Modify files in `Dataset_Structural_GPCRs` that belong to the published v1.0.0 snapshot without
  an explicit request. Data changes go to a new version (v1.1.0).
- Commit large data or run outputs to this repository.
- Push, merge to `main`, or open PRs unless asked.
- Change objective definitions, operator behavior, or default hyperparameters silently: they
  determine published results. Propose the change first.
- Modify the irace tuning code (`org.tm_msaligner.auto`) or the legacy non-structural classes it
  depends on (`TM_MSASolution`, `SPXMSACrossover`, the old mutations, `MultiObjTMMSAProblem`,
  etc.). The authors will decide later whether to continue that line. Only API renames that a
  dependency upgrade makes unavoidable are allowed.

## 9. Future work: graphical user interface

Once the rewrite is complete, the authors plan a GUI for the tool, probably a Python
**Streamlit** application that drives the Java core. Design the core so that a GUI can be added
without touching the algorithm:

- Keep all logic in the library (`org.tmm2structalign`); `cli` is a thin layer. No
  `System.exit`, no console-only output, and no interactive prompts in library code.
- Every run must be fully described by a configuration file (`.properties` or JSON) and must
  write **machine-readable** outputs with a stable layout: `config.json`, `metadata.json`,
  `FUN.tsv`, `MSASol*.fasta`, and a progress log that can be followed while the run is active
  (one JSON object per line, e.g. `progress.jsonl` with evaluations and best objective values).
- Report errors with specific exceptions and clear messages that the CLI turns into a non-zero
  exit code plus a JSON error file, so a GUI can display them.
- Keep the CLI arguments and output formats stable and documented; treat them as the public API
  between the Java core and the GUI.
- Do not start the GUI before the rewrite (Part 1 of the plan) is finished.

