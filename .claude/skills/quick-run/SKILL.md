---
name: quick-run
description: Build the Transept core jar and run a short alignment of classF (a few seconds) to check a change end to end, then summarise the front and the files of the run. Use after changing the search, the objectives, the operators or the outputs, or when asked to try the tool.
---

# A quick run of Transept core

1. Find the benchmark data. Transept core expects `../Dataset_Structural_GPCRs/GPCRdb`. If it is
   missing, tell the user how to get it and stop:

   ```bash
   git clone https://github.com/TM-M2StructAlign/Dataset_Structural_GPCRs.git ../Dataset_Structural_GPCRs
   git -C ../Dataset_Structural_GPCRs checkout v1.0.0
   ```

2. Build the jar without tests (the `check` skill runs them):

   ```bash
   mvn -q clean package -DskipTests
   ```

3. Run classF (11 proteins) with a small budget, in a fresh directory outside the repository:

   ```bash
   OUT=$(mktemp -d)
   java -jar target/transept-core-0.1.0-jar-with-dependencies.jar run \
     --data-root ../Dataset_Structural_GPCRs/GPCRdb --dataset classF \
     --seeds clustalw,kalign,mafft,muscle5,tcoffee --repair-seeds \
     --output-root "$OUT" --run-id "quick_$(date +%Y%m%d_%H%M%S)" --seed 1 \
     --max-evaluations 500 --population-size 20
   ```

   Add the options being tested (e.g. `--algorithm`, `--mutation`, `--objectives`). Use
   `--seeds internal` to try a dataset without precomputed alignments.
4. Summarise the run directory `$OUT/ejecuciones/classF/<run-id>/`:
   - if `error.json` exists, show its message and stop;
   - `FUN.csv` (no header, one row per solution, objectives in the order of `--objectives`): the
     number of solutions and the best value of each objective;
   - `runtime.csv`: the time taken;
   - `metadata.json`: the build, the seeds used and the mapping status of the structures (report
     any structure that was not mapped).
5. Leave the directory in the temporary folder and say where it is; do not copy run outputs into
   the repository.
