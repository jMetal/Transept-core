# Output files

## Run outputs

Every `run` writes to `<output-root>/ejecuciones/<dataset>/<run-id>/` and refuses to start if the
directory exists. The layout and formats are stable: they are the interface for scripts and for
the planned graphical interface.

| File | Content |
| --- | --- |
| `FUN.csv` | One line per solution of the result of the algorithm: the objectives (by default SOPwTP and lDDT; see `--objectives`), comma-separated, no header |
| `MSASol<i>.fasta` | Alignment of line `i` of `FUN.csv` (0-based), aligned FASTA |
| `config.json` | Full configuration of the run |
| `metadata.json` | Provenance and diagnostics of the run |
| `progress.jsonl` | One JSON object per line, every `--progress-interval` evaluations |
| `runtime.csv` | Header `runtime_ms,runtime_seconds,evaluations` and one row |
| `seed.txt` | Random seed |
| `seeds/` | The internal seeds computed by the run, if it computed them (not when it read them from `internal_seeds/`) |
| `error.json` | Only if the run failed: exception type and message |

### `FUN.csv`

Both objectives are written as computed, and **higher is better for both**. jMetal minimises, so
the solutions store them negated internally, but the file does not.

The solutions are the result returned by the algorithm: the final population (NSGA-II, NSGA-III,
AGE-MOEA, SMS-EMOA, MOEA/D, the genetic algorithm, random search) or the archive (PAES, and
NSGA-II configured with `--algorithmResult externalArchive`). The result is not filtered, so
it can contain dominated solutions and repeated ones (two lines with the same values, usually the
same alignment); `scripts/pareto_analysis.py` keeps the non-dominated ones (repeated lines included). Example:

```text
61113.0,0.565144899719245
48650.0,0.6140871258624288
```

### `config.json`

Every option of the run, with defaults resolved, for example
`"mutationType":"TM_SAFE_GAP_SHIFT"`, `"mutationProbability":0.2` and `"lddtCutoff":15.0`.
`configuration` is the configuration string given by the user (empty if none), and the stopping
condition is `maxEvaluations` or `maxTimeSeconds` (the other one is `null`). Together with
the dataset version and the git commit of `metadata.json`, it is enough to repeat the run.

### `metadata.json`

| Key | Content |
| --- | --- |
| `tool`, `algorithmType`, `algorithmName` | What ran |
| `configuration` | For algorithms built from a parameter space (NSGA-II), the complete configuration used, `--parameter value ...`; empty for the others |
| `build` | Git commit and branch, whether the working tree had uncommitted changes, Java version |
| `startedAt`, `finishedAt`, `runtimeMs`, `evaluations`, `frontSize` | Timing and size of the result |
| `stoppedByUser` | `true` if the run was stopped early with a `STOP` file ([cli](cli.md#stopping-a-run-early)): `evaluations` is then below the budget, and the run cannot be repeated exactly |
| `objectives` | Objective names, in the order of `FUN.csv` |
| `seeds` | For every seed: `method` (`internal:biojava-g<open>x<extension>` for internal seeds), file, SHA-256, `duplicateOf` if it is identical to an earlier seed, `repairedGapRuns` (gap runs moved out of helices by the repair) and `origin` (`file`, `precomputed` or `computed`) |
| `structures` | For every protein: mapping `status`, sequence `identity` and `coverage` of the mapping, and a `detail` message |
| `proteinsWithoutStructure` | Proteins whose structure is not used by the lDDT |

Mapping statuses: `NUMBERED` (PDB numbering, residues checked), `ALIGNED` (numbering disagreed;
mapped by sequence alignment), `UNVERIFIED` (numbering used without checks), and three statuses
whose structure is not used: `LOW_IDENTITY`, `MISSING_FILES` and `INVALID_FILES`.

### `progress.jsonl`

```json
{"evaluations":200,"elapsedMs":721,"best":{"SOPwTP":60715.0,"lDDT":0.5897227517717524}}
```

`best` holds the best value of each objective in the current result of the search (population or
archive), each possibly from a
different solution. The file can be followed while the run is active.

## `baseline` output

CSV with the header `dataset,method,SOPwTP,lDDT,file,status,message`, one row per precomputed
alignment. `status` is `ok` or `failed`; failed rows leave the objectives empty and explain
the problem in `message`.

## Evaluation output

`score` writes one row per alignment with the columns `relative_path`, `run_id`, `msa_file`,
`sequence_count`, the precision, recall and F1 of all pairs and of TM pairs, the TM-gap rate, the raw
counts behind them (`reference_pairs`, `test_pairs`, `correct_pairs`, the same for TM pairs,
`tm_gap_count`, `tm_gap_opportunities`), `status` and `message`. `status` is `ok` or `failed`; failed
rows leave the metrics empty and explain the problem in `message`. See [metrics](metrics.md).
