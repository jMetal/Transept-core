# AGENTS.md

Guidance for AI coding assistants (and humans) working on **TM-M2StructAlign**. Read this file
first, then `JAVA_CODING_GUIDELINES.md` and `GIT_GUIDELINES.md`.

## 1. What this project is

TM-M2StructAlign is a multiobjective evolutionary algorithm (NSGA-II built from jMetal
components) for multiple sequence alignment (MSA) of transmembrane proteins, specifically human
GPCRs. It optimizes two objectives:

1. **SOPwTP** (`objective.TopologyAwareSumOfPairs`): sum of pairs using PHAT when both residues
   are transmembrane (TM) and BLOSUM62 otherwise, with affine gap penalties that differ for TM
   (8/3) and non-TM (3/1).
2. **lDDT** (`objective.StructuralLddt`): reference-free, superposition-free lDDT computed from
   AlphaFold2 Cα distance matrices (R0 = 15 Å).

Both are maximized; `problem.StructuralMsaProblem` negates them for jMetal minimization.
`--objectives` selects other formulations for experiments (the lDDT alone, or a third objective,
the aligned-pair coverage); the code does not assume two objectives, except the Python scripts
and the studio, which expect the default.

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
mvn clean verify     # unit + integration tests, Checkstyle, JaCoCo; must pass before committing
mvn clean package    # target/TM-MSAStructAligner-1.0-jar-with-dependencies.jar
GPCR_DATASET=../Dataset_Structural_GPCRs/GPCRdb mvn verify   # also runs real-data ITs
```

The JAR runs the picocli command line (`cli.TmStructAlignCli`):

```bash
java -jar target/TM-MSAStructAligner-1.0-jar-with-dependencies.jar run \
  --data-root ../Dataset_Structural_GPCRs/GPCRdb --dataset classF \
  --seeds clustalw,kalign,mafft,muscle5,tcoffee,tmaligner \
  --output-root results --run-id seed_000 --seed 0 \
  --max-evaluations 10000 --population-size 100 --threads 4
java -jar target/TM-MSAStructAligner-1.0-jar-with-dependencies.jar baseline \
  --data-root ../Dataset_Structural_GPCRs/GPCRdb --datasets classF --output baselines.csv
```

Run outputs (`<output-root>/ejecuciones/<dataset>/<run-id>/`): `FUN.csv` (no header),
`MSASol<i>.fasta`, `runtime.csv`, `config.json`, `metadata.json`, `progress.jsonl`, `seed.txt`,
and `error.json` on failure. Existing run directories are never overwritten. Tabular outputs are
CSV; `scripts/pareto_analysis.py` also reads `FUN.tsv`/`runtime.txt` of the original runner.

For quick experiments use `classF` (11 sequences) with ~500 evaluations and a population of 20:
it runs in a few seconds. Tests use a generated dataset (`testing.SyntheticDataset`).

## 4. Code map (`src/main/java/org/tmm2structalign/`)

| Package | Role |
| --- | --- |
| `model` | `ProteinSequence` (residues + topology), `ProteinSet`, `Topology` |
| `alignment` | `GapIntervals` (immutable gap runs of one row), `Msa` (decoded integer view, TM flag of gaps), `GapColumnCompactor`, `MsaSolution` (jMetal solution, one `GapIntervals` per protein) |
| `io` | Readers for 3-line topology, aligned FASTA (tool-specific ids), `_D.csv`/`_idx.csv`, PDB chains; FASTA writer |
| `structure` | `SequenceStructureMapper` (numbering checked by residue identity, alignment fallback), `StructureLoader`, `ProteinStructures`, `StructuralNeighborhood` |
| `objective` | `TopologyAwareSumOfPairs`, `StructuralLddt`, `AlignedPairCoverage`, `ObjectiveSet` (formulations of `--objectives`), `SubstitutionMatrix` |
| `operator` | `TmSafeSinglePointCrossover`, mutations that never break a helix (`TmSafeGapShiftMutation` (default), `MergeAdjacentGapRunsMutation`, `RandomChoiceMutation`), `TmBoundaries`, `RandomSource`. Operators evaluated and removed are listed in tutorial T2 and the internal notes |
| `problem` | `StructuralMsaProblem` |
| `algorithm` | `SeededPopulationCreation`, `AlgorithmConfiguration`, `AlgorithmType`, `AlgorithmFactory` (NSGA-II, SPEA2 and the GA from `evolver`; asynchronous NSGA-II and GA from `jmetal-parallel`; NSGA-III, AGE-MOEA, SMS-EMOA, MOEA/D, PAES and random search built with jMetal), `SearchRun`, `RandomAlignmentCreation`, `PerturbedSeedCreation` |
| `evolver` | Algorithms built with Evolver from a parameter space (`src/main/resources/parameterSpaces/*.yaml`) and a default configuration (`defaultConfigurations/*Default.txt`): `MsaNsga2`, `MsaSpea2`, `MsaGeneticAlgorithm`, their parameters (`MsaParameterFactory`, `Msa*Parameter`), `MsaComponents`, `MsaOperators`, `MsaConfigurations`, `MsaRunContext` |
| `run` | `RunConfiguration`, `Experiment`, `ExperimentRunner`, `ProgressObserver`, `BaselineScorer`, `Json` |
| `cli` | `run` and `baseline` commands |
| `evaluation` | Comparison with GPCRdb reference alignments (`EvaluationMain`) |

`org.tm_msaligner` contains only the irace tuning code (`auto`) and the 27 legacy classes it
needs; see §9.

## 5. Conventions

- Follow `JAVA_CODING_GUIDELINES.md`. Platform: Java 21, jMetal 7.10, JUnit 6 (Jupiter) +
  Mockito (`mock`/`when`/`verify` only), unit tests `*Test` (surefire) and integration tests
  `*IT` (failsafe, `mvn verify`).
- Checkstyle (Google style) and JaCoCo (>= 80 % line coverage) are enforced on
  `org.tmm2structalign` only; `config/checkstyle/suppressions.xml` allows the `IT` suffix.
- Branch `feat/algorithms` builds against **jMetal 7.11-SNAPSHOT** (the SMS-EMOA, MOEA/D and PAES
  builders accept a given initial population or solution only from 7.11). Install the snapshot
  from the local jMetal sources until 7.11 is on Maven Central, then switch to `7.11` and merge.
  The asynchronous algorithms also need the 7.11-SNAPSHOT of `jmetal-parallel` (initial
  population of the asynchronous genetic algorithms).
- NSGA-II, SPEA2 and the GA are built with **Evolver** (`org.uma.jmetal:Evolver:3.0-SNAPSHOT`, GPL-3): install it
  with `mvn install -DskipTests` from the `develop` branch of the sibling `Evolver` repository.
  Their parameter spaces (`MSA*.yaml`) and default configurations (`MSA*Default.txt`) are
  part of the public API (Studio reads them from the JAR). The licence of the combination is
  pending a decision of the authors.
- Follow `GIT_GUIDELINES.md`: Conventional Commits, one logical change per commit, never mix
  `test:` with `fix:`/`feat:`, never mix code with docs; every commit must build.
- Write code, comments, and messages in English. Discussion and review documents may be in
  Spanish.
- Draw all randomness from `operator.RandomSource` (default: the global `JMetalRandom`) so that
  runs are reproducible from `--seed`; tests inject `testing.SeededRandom` or Mockito mocks.

## 6. Pitfalls (read before touching these areas)

1. **Objective 2 differs from the published one.** The lDDT is the mean of the anchor scores; the
   original aggregation (harmonic mean with occupancy, divided by the length) rewarded compact
   alignments and was removed (finding A3). Like any superposition-free local score, it is nearly
   blind to helix register shifts: check alignments against reference-based metrics.
2. **Structure mapping.** The verified mapping corrects OXER1 and MCHR1, which the original
   implementation mis-mapped; it changes the lDDT of classA_002, classA_004 and classA_004_19.
   The published values can only be reproduced with tag `legacy-v1.0`.
3. **Seeds.** Caretta outputs are residue-level misaligned for most classes, and the structural
   tools' outputs in the full classA_002/003/004 datasets do not match the GPCRdb sequences.
   Choose `--seeds` deliberately.
4. **Differential oracle.** `BaselineOracleIT` checks SOPwTP against the original values of the
   dataset's `baselines.csv` and lDDT against a regression snapshot (`lddt-snapshot.csv`), both in
   `src/test/resources/oracle`. If you change an objective on purpose, regenerate the snapshot in
   the same change and explain why.
5. `.gitignore` excludes `*.fasta` except under `src/test/resources/`. Prefer generating fixtures
   in `@TempDir`.
6. **TM integrity is not a formal constraint** (`numberOfConstraints() = 0`), but no operator
   breaks a helix: the crossover never cuts one and every mutation forbids it. Gaps inside helices
   inherited from the seeds are kept. Do not change defaults without asking the authors.
7. **Superseded components are removed, not kept as options** (the original mutation, the
   original lDDT aggregation, the legacy mapping, the two-point crossover, the original filling
   of the initial population): reproducing the published results is the job of tag
   `legacy-v1.0`.
8. Results are not bit-identical to the published ones (jMetal 6.2 → 7.10 changes how random
   numbers are consumed); same seed + same version is bit-identical.

## 7. Current work (November 2026)

Branch `antonio` (and `feat/algorithms` for the algorithms and operators). The authors' internal
review notes (code review with ids A*/M*/B*, dataset review with ids D*, work plan, operator and
initialization proposals, diagnostic scripts) are kept **outside the repository**, in the sibling
directory `../notes/`, and are not versioned; versions up to 2026-11-04 remain in the git history
under `docs/review/`. Do not add internal notes to this repository.
Part 1 of the plan (platform migration, rewrite of the structural core with tests, removal of
the old structural flow) is done; tag `legacy-v1.0` marks the pre-rewrite state. Part 2
(new operators, objectives, algorithm improvements) is pending; Phase 0 decisions of the plan
belong to the authors, so ask before implementing work that depends on them.

When fixing a finding, reference its id in the commit body (e.g. `Refs: A2`) and, if the notes
are available locally, update `../notes/review/02-code-review.md`.

## 8. Documentation

User documentation is Markdown in the repository (no site generator): the short `README.md`,
tutorials in `docs/tutorials/` (T1, T2, ...) and reference pages in `docs/reference/` (`cli.md`,
`outputs.md`, `metrics.md`). Write it in English. Each tutorial starts with a table (level,
version, time, machine of the timings, prerequisites) and has numbered steps with the exact
commands and their expected output.

- When you add or change an option, an output file or an operator, update `docs/reference/` and
  the affected tutorials in a `docs:` commit.
- `example.FormulationTour` prints the example of T2; `FormulationTourTest` checks that T2 shows
  its current output. If the output changes, paste the new output into T2.
- Internal notes do not belong here (see §7).

## 9. Do not

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

## 10. Graphical user interface

A first prototype exists in the sibling repository `TM-M2StructAlign-Studio` (Streamlit, local
use). It drives this tool only through the `run` command and the output files, so treat the CLI
options and the output formats as its public API: when you change them, update
`tmstudio/options.py` and the readers of the studio.

The original guidance for the interface follows.

Once the rewrite is complete, the authors plan a GUI for the tool, probably a Python
**Streamlit** application that drives the Java core. Design the core so that a GUI can be added
without touching the algorithm:

- Keep all logic in the library (`org.tmm2structalign`); `cli` is a thin layer. No
  `System.exit`, no console-only output, and no interactive prompts in library code.
- Every run must be fully described by a configuration file (`.properties` or JSON) and must
  write **machine-readable** outputs with a stable layout: `config.json`, `metadata.json`,
  `FUN.csv`, `MSASol*.fasta`, and a progress log that can be followed while the run is active
  (one JSON object per line, e.g. `progress.jsonl` with evaluations and best objective values).
- Report errors with specific exceptions and clear messages that the CLI turns into a non-zero
  exit code plus a JSON error file, so a GUI can display them.
- Keep the CLI arguments and output formats stable and documented; treat them as the public API
  between the Java core and the GUI.
- Do not start the GUI before the rewrite (Part 1 of the plan) is finished.

