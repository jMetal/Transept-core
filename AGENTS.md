# AGENTS.md

Guidance for AI coding assistants (and humans) working on **Transept core**, the engine of
[Transept](https://github.com/jMetal/Transept). Read this file
first, then `JAVA_CODING_GUIDELINES.md` and `GIT_GUIDELINES.md`.

## 1. What this project is

Transept core is a multiobjective evolutionary tool (NSGA-II by default, built with jMetal and
Evolver) for multiple sequence alignment (MSA) of transmembrane proteins, specifically human
GPCRs. It derives from the TM-M2StructAlign software. By default it optimizes two objectives:

1. **SOPwTP** (`objective.TopologyAwareSumOfPairs`): sum of pairs using PHAT when both residues
   are transmembrane (TM) and BLOSUM62 otherwise, with affine gap penalties that differ for TM
   (8/3) and non-TM (3/1).
2. **lDDT** (`objective.StructuralLddt`): reference-free, superposition-free lDDT computed from
   AlphaFold2 Cα distance matrices (R0 = 15 Å).

Both are maximized; `problem.StructuralMsaProblem` negates them for jMetal minimization.
`--objectives` selects other formulations for experiments (the lDDT alone, or a third objective,
the aligned-pair coverage); the code does not assume two objectives.

## 2. Workspace layout

The benchmark data is expected in a sibling directory:

```text
<parent>/
├── Transept-core/               # this repository (Java, Maven)
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
mvn clean package    # target/transept-core-0.1.0-jar-with-dependencies.jar
GPCR_DATASET=../Dataset_Structural_GPCRs/GPCRdb mvn verify   # also runs real-data ITs
```

The JAR runs the picocli command line (`cli.TranseptCli`):

```bash
java -jar target/transept-core-0.1.0-jar-with-dependencies.jar run \
  --data-root ../Dataset_Structural_GPCRs/GPCRdb --dataset classF \
  --seeds clustalw,kalign,mafft,muscle5,tcoffee,tmaligner \
  --output-root results --run-id seed_000 --seed 0 \
  --max-evaluations 10000 --population-size 100 --threads 4
java -jar target/transept-core-0.1.0-jar-with-dependencies.jar baseline \
  --data-root ../Dataset_Structural_GPCRs/GPCRdb --datasets classF --output baselines.csv
```

Run outputs (`<output-root>/ejecuciones/<dataset>/<run-id>/`): `FUN.csv` (no header),
`MSASol<i>.fasta`, `runtime.csv`, `config.json`, `metadata.json`, `progress.jsonl`, `seed.txt`,
and `error.json` on failure. Existing run directories are never overwritten. Tabular outputs are
CSV; `scripts/pareto_analysis.py` also reads `FUN.tsv`/`runtime.txt` of the original runner.

Campaigns on a Slurm cluster: `scripts/slurm/prepare_campaign.sh <template> <output>` builds a
self-contained package (jar, provenance, data subset, scripts) from a template such as
`scripts/slurm/formulations/` (`cells.tsv`, `configurations.tsv`, `campaign.env`, `README.md`);
its `submit.sh`, `check.sh` and `analyse.sh` run it unattended and resumably. Keep the packages
outside the repository (e.g. `../campaigns/`).

For quick experiments use `classF` (11 sequences) with ~500 evaluations and a population of 20:
it runs in a few seconds. Tests use a generated dataset (`testing.SyntheticDataset`).

## 4. Code map (`src/main/java/org/transept/`)

| Package | Role |
| --- | --- |
| `model` | `ProteinSequence` (residues + topology), `ProteinSet`, `Topology` |
| `alignment` | `GapIntervals` (immutable gap runs of one row), `Msa` (decoded integer view, TM flag of gaps), `GapColumnCompactor`, `MsaSolution` (jMetal solution, one `GapIntervals` per protein) |
| `io` | Readers for 3-line topology, aligned FASTA (tool-specific ids), `_D.csv`/`_idx.csv`, PDB chains; FASTA writer |
| `structure` | `SequenceStructureMapper` (numbering checked by residue identity, alignment fallback), `StructureLoader`, `ProteinStructures`, `StructuralNeighborhood` |
| `objective` | `TopologyAwareSumOfPairs`, `StructuralLddt`, `AlignedPairCoverage`, `ObjectiveSet` (formulations of `--objectives`), `SubstitutionMatrix` |
| `operator` | `TmSafeSinglePointCrossover`, mutations that never break a helix (`TmSafeGapShiftMutation` (default), `MergeAdjacentGapRunsMutation`, `RandomChoiceMutation`), `TmHelixRepair`, `TmBoundaries`, `RandomSource`. Operators evaluated and removed are listed in tutorial T2 |
| `problem` | `StructuralMsaProblem` |
| `algorithm` | `SeededPopulationCreation`, `AlgorithmConfiguration`, `AlgorithmType`, `AlgorithmFactory` (NSGA-II, SPEA2 and the GA from `evolver`; asynchronous NSGA-II and GA from `jmetal-parallel`; NSGA-III, AGE-MOEA, SMS-EMOA, MOEA/D, PAES and random search built with jMetal), `SearchRun`, `RandomAlignmentCreation`, `PerturbedSeedCreation` |
| `evolver` | Algorithms built with Evolver from a parameter space (`src/main/resources/parameterSpaces/*.yaml`) and a default configuration (`defaultConfigurations/*Default.txt`): `MsaNsga2`, `MsaSpea2`, `MsaGeneticAlgorithm`, their parameters (`MsaParameterFactory`, `Msa*Parameter`), `MsaComponents`, `MsaOperators`, `MsaConfigurations`, `MsaRunContext` |
| `run` | `RunConfiguration`, `Experiment`, `ExperimentRunner`, `ProgressObserver`, `BaselineScorer`, `Json` |
| `cli` | `run` and `baseline` commands |
| `evaluation` | Comparison with GPCRdb reference alignments (`EvaluationMain`) |

## 5. Conventions

- Follow `JAVA_CODING_GUIDELINES.md`. Platform: Java 21, jMetal 7.11, JUnit 6 (Jupiter) +
  Mockito (`mock`/`when`/`verify` only), unit tests `*Test` (surefire) and integration tests
  `*IT` (failsafe, `mvn verify`).
- Checkstyle (Google style) and JaCoCo (>= 80 % line coverage) are enforced on
  `org.transept`; `config/checkstyle/suppressions.xml` allows the `IT` suffix.
- The build uses **jMetal 7.11**, from Maven Central: it is the first version whose SMS-EMOA,
  MOEA/D and PAES builders and asynchronous genetic algorithms (`jmetal-parallel`) accept a given
  initial population or solution.
- NSGA-II, SPEA2 and the GA are built with **Evolver 2.6** (`org.uma.jmetal:Evolver`, GPL-3,
  from Maven Central). Their parameter spaces (`MSA*.yaml`) and default configurations (`MSA*Default.txt`) are
  part of the public API (the Transept interface reads them from the JAR).
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
   The published values can only be reproduced with the original TM-M2StructAlign software
   (tag `legacy-v1.0` of its repository).
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
   of the initial population): reproducing the published results is the job of the original
   TM-M2StructAlign software.
8. Results are not bit-identical to the published ones (jMetal 6.2 → 7.10 changes how random
   numbers are consumed); same seed + same version is bit-identical.

## 7. Documentation

User documentation is Markdown in the repository (no site generator): the short `README.md`,
tutorials in `docs/tutorials/` (T1, T2, ...) and reference pages in `docs/reference/` (`cli.md`,
`outputs.md`, `metrics.md`). Write it in English. Each tutorial starts with a table (level,
version, time, machine of the timings, prerequisites) and has numbered steps with the exact
commands and their expected output.

- When you add or change an option, an output file or an operator, update `docs/reference/` and
  the affected tutorials in a `docs:` commit.
- `example.FormulationTour` prints the example of T2; `FormulationTourTest` checks that T2 shows
  its current output. If the output changes, paste the new output into T2.
- Internal notes and review documents do not belong in this repository.

## 8. Do not

- Modify files in `Dataset_Structural_GPCRs` that belong to the published v1.0.0 snapshot without
  an explicit request. Data changes go to a new version (v1.1.0).
- Commit large data or run outputs to this repository.
- Push, merge `develop` into `main`, or open PRs unless asked (the working branch is `develop`).
- Change objective definitions, operator behavior, or default hyperparameters silently: they
  determine published results. Propose the change first.
## 9. Graphical user interface

The [Transept](https://github.com/jMetal/Transept) interface (Streamlit) drives this tool only
through its commands (`run`, `seeds`, `algorithms`) and the output files, and reads the parameter
spaces and default configurations from the JAR. Treat the CLI options, the output formats and
those resources as its public API: when you change them, update the interface in step.

- Keep all logic in the library (`org.transept`); `cli` is a thin layer. No `System.exit`, no
  console-only output, and no interactive prompts in library code.
- Every run is fully described by its configuration and writes machine-readable outputs with a
  stable layout (`config.json`, `metadata.json`, `FUN.csv`, `MSASol*.fasta`, `progress.jsonl`).
- Report errors with specific exceptions and clear messages that the CLI turns into a non-zero
  exit code plus a JSON error file, so the interface can display them.
