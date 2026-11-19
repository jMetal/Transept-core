---
name: cli-contract
description: Keep Transept core and the Transept interface in step when a command, an option, an option value, an output file or a parameter space of the core changes. Use whenever such a change is planned or made in Transept core.
---

# Change the contract between Transept core and Transept

The Transept interface (repository `jMetal/Transept`) drives this tool only through its commands,
its output files and the parameter spaces and default configurations inside the jar. A change to
any of them is a change of contract.

## In Transept core (this repository)

1. Make the code change and its tests (separate `feat:`/`fix:` and `test:` commits).
2. Update the documentation in a `docs:` commit:
   - `docs/reference/cli.md` for commands, options and their values or defaults;
   - `docs/reference/outputs.md` for files of a run (`FUN.csv`, `config.json`, `metadata.json`,
     `progress.jsonl`, `error.json`, ...) and reports of `baseline` and `score`;
   - `docs/reference/configuration.md` and `docs/reference/operators.md` for parameter spaces
     (`src/main/resources/parameterSpaces/`) and operators;
   - the tutorials of `docs/tutorials/` that show the changed option or output. If the output of
     `example.FormulationTour` changes, paste it into T2 (`FormulationTourTest` checks it).
3. Add the change to `CHANGELOG.md` under `[Unreleased]`.
4. Run the `check` skill.

## In Transept (the interface)

Look for a clone next to this one (`../Transept` or `../transept`). If there is none, list what must
change there and stop.

1. Options and values: `transept/options.py`, and the pages that offer them (`pages/new_run.py`).
2. Readers of the outputs: `transept/runs.py`, `runner.py`, `reference_scores.py`, `baselines.py`,
   `internal_seeds.py`, `repair.py`, `algorithms.py`, `spaces.py`.
3. Tests and fixtures: `tests/conftest.py` (the toy runs must look like real ones) and the tests of
   the modules changed.
4. Rebuild the jar of the interface with `make core` and run `make lint` and `make test`;
   `tests/test_options_drift.py` compares `options.py` with the help of the jar and must pass, not
   be skipped.
5. Older runs must still open: if an output changes, keep reading the previous layout.

Each repository gets its own commits and its own pull request into `develop`; link the two pull
requests to each other. Do not push or open them unless the user asks.
