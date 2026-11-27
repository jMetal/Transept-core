# T3. A Study With Several Seeds

| | |
| --- | --- |
| **Level** | Intermediate |
| **Version** | 1.0 (2026-11-04) |
| **Time** | about 20 minutes, plus 10 runs of about 10 seconds each |
| **Timings measured on** | Apple M5 Pro (18 cores), 64 GB of RAM, macOS 26.6.2, Java 21.0.12 (Oracle JDK), Python 3.9 with NumPy and Matplotlib 3.9 |
| **Prerequisites** | [T1. First run in five minutes](T1-first-run.md), [T4. Choosing an alignment](T4-choosing-and-evaluating.md) |

The search is stochastic: two runs with different random seeds return different alignments. A
study therefore runs the same configuration with several seeds and analyses the runs together.
This tutorial runs 10 seeds on `classF` with a reduced budget, combines their fronts, measures how
much the runs differ, and picks alignments from the combined front.

## Step 1: run 10 seeds

On Linux or macOS:

```bash
for seed in 0 21 42 63 84 105 126 147 168 189; do
  java -jar target/transept-core-1.0-SNAPSHOT-jar-with-dependencies.jar run \
    --data-root ../Dataset_Structural_GPCRs/GPCRdb --dataset classF \
    --seeds clustalw,kalign,mafft,muscle5,tcoffee,tmaligner \
    --output-root results --run-id seed_$(printf %03d $seed) --seed $seed \
    --max-evaluations 2500 --population-size 100 --threads 4
done
```

On Windows, `scripts/run_10_seeds_windows.ps1` does the same with the same seeds (run identifiers
`run_01` ... `run_10`), and checks that every run wrote its files:

```powershell
.\scripts\run_10_seeds_windows.ps1 -JarPath target\transept-core-1.0-SNAPSHOT-jar-with-dependencies.jar `
  -DataDirectory ..\Dataset_Structural_GPCRs\GPCRdb -ProblemName classF -OutputRoot results `
  -SeedMethods clustalw,kalign,mafft,muscle5,tcoffee,tmaligner -MaxEvaluations 2500
```

Each run takes about 10 seconds. The budget of 2 500 evaluations keeps the tutorial short; a real
study uses more (Step 6).

## Step 2: analyse the runs together

`scripts/pareto_analysis.py` reads every run of a dataset under `results/ejecuciones/` (it needs
Python 3 with NumPy and Matplotlib):

```bash
python3 scripts/pareto_analysis.py results classF --output results/analysis_classF
```

```text
Analyzed 1000 solutions from 10 runs
Combined non-dominated objective vectors: 67
Balanced compromise: seed_063/MSASol23.fasta
Results: results/analysis_classF
```

It writes:

| File | Content |
| --- | --- |
| `combined_solutions.csv` | Every solution of every run, with its objectives normalised over all the runs |
| `combined_pareto_front.csv` | The non-dominated solutions of the union of the runs |
| `indicators_by_run.csv` | Quality indicators of every run (Step 3) |
| `indicator_summary.csv` | Median and interquartile range of each indicator over the runs |
| `selected_solutions.csv`, `selected_alignments/` | The alignments chosen from the combined front (Step 4) |
| `pareto_front_classF.{png,svg,pdf}` | The figure below |

![Combined front of 10 runs on classF](figures/T3-pareto-front-classF.png)

The grey points are the 1 000 final solutions; the blue line is the combined non-dominated front.
Both objectives are scaled to [0, 1] over all the runs, so that they can be compared.

## Step 3: how much do the runs differ?

`indicators_by_run.csv` measures every run against the combined front:

| Indicator | Meaning | Better |
| --- | --- | --- |
| `hypervolume` | Area dominated by the run's front in the normalised objective space | Higher |
| `igd_plus` | Distance from the combined front to the run's front | Lower |
| `nondominated_solutions` | Distinct non-dominated solutions of the run | |
| `combined_front_contribution` | Solutions of the run that are on the combined front | |
| `reached_compromise_region` | Whether the run has a solution almost as close to the ideal point as the balanced compromise (within `--epsilon`, 0.05) | |

```text
run       hypervolume  igd_plus  nondominated  contribution
seed_000  0.551        0.039     45             10
seed_021  0.464        0.093     43             0
seed_042  0.499        0.072     50             3
seed_063  0.554        0.031     41             48
seed_084  0.500        0.075     61             0
seed_105  0.533        0.051     59             0
seed_126  0.512        0.068     62             0
seed_147  0.471        0.084     48             1
seed_168  0.550        0.041     53             6
seed_189  0.609        0.009     55             64
```

The hypervolume goes from 0.46 to 0.61 (median 0.52, interquartile range 0.05). Two runs
(`seed_063` and `seed_189`) provide most of the combined front, and four contribute nothing to
it. A conclusion drawn from a single run
would depend heavily on the seed: this is why studies use 10 or more, and report medians and
spreads instead of single values.

## Step 4: choose alignments from the combined front

`selected_solutions.csv` lists four alignments of the combined front, copied to
`selected_alignments/`:

| Selection | Rule | Run / row | SOPwTP | lDDT |
| --- | --- | --- | --- | --- |
| `extreme_f1` | Best SOPwTP | seed_147 / 0 | 61 372 | 0.8004 |
| `extreme_f2` | Best lDDT | seed_189 / 1 | 55 407 | 0.8234 |
| `compromise_balanced` | Closest to the ideal point with weights (0.5, 0.5) | seed_063 / 23 | 58 856 | 0.8152 |
| `compromise_sensitivity` | The same with weights (0.7, 0.3) | seed_189 / 48 | 60 261 | 0.8088 |

The weights are options of the script (`--weight-f1`, `--weight-f2`, and
`--sensitivity-weight-f1`, `--sensitivity-weight-f2`). Giving more weight to SOPwTP moves the
compromise towards the best-SOPwTP end of the front.

Evaluate them against the GPCRdb reference as in [T4](T4-choosing-and-evaluating.md):

```bash
java -jar target/transept-core-1.0-SNAPSHOT-jar-with-dependencies.jar score \
  --reference ../Dataset_Structural_GPCRs/GPCRdb/reference_alignments/classF.fasta \
  --topology ../Dataset_Structural_GPCRs/GPCRdb/sequences/tmregions/classF_predicted_topologies.3line \
  --msa-root results/analysis_classF/selected_alignments --glob '*.fasta' \
  --output results/evaluation_selected_classF.csv
```

```text
compromise_balanced.fasta     Pair-F1=0.556866  TM-Pair-F1=0.999820  TM-gap-rate=0.005170
compromise_sensitivity.fasta  Pair-F1=0.554828  TM-Pair-F1=0.999820  TM-gap-rate=0.005170
extreme_f1.fasta              Pair-F1=0.540159  TM-Pair-F1=0.999820  TM-gap-rate=0.005170
extreme_f2.fasta              Pair-F1=0.568437  TM-Pair-F1=0.999820  TM-gap-rate=0.004318
```

As in T4, the agreement with the reference grows from the best-SOPwTP end to the best-lDDT end:
the best-lDDT extreme is the best of the four, and the best-SOPwTP extreme the worst.

## Step 5: compare with the seeds

The combined front is only interesting if it improves on the starting alignments. Evaluate the
seeds as in [T4, Step 5](T4-choosing-and-evaluating.md#step-5-compare-with-the-seeds): on
`classF` the balanced compromise (pair F1 0.5569) is better than five of the six seeds and
slightly worse than ClustalW (0.5598), and the best-lDDT extreme (0.5684) is better than all of
them.

## Step 6: towards a real study

This tutorial uses a small budget and one small dataset. A study meant for publication should
also:

- **Use a larger budget** (`--max-evaluations`), and check with `progress.jsonl` that the best values
  have stopped improving before the end.
- **Run more datasets.** Runtime grows with the square of the number of proteins and with their
  length: the `_19` variants (19 proteins) are a good intermediate step before the full classes.
- **Fix everything in advance and keep the metadata.** Every run records its configuration
  (`config.json`) and provenance (`metadata.json`: git commit, SHA-256 of the seeds, mapping of
  the structures). Keep them with the results.
- **Compare against baselines**, such as the seeds themselves and the random searches of
  [T5](T5-comparing-algorithms-and-operators.md).

This version does not reproduce the published experiments: it corrects the structure mapping of
OXER1 and MCHR1, removes gap-only columns, defines the structural objective as the mean of the
anchor scores, uses a mutation that never breaks a helix and fills the initial population without
duplicates. To reproduce the published results,
use the original TM-M2StructAlign at tag `legacy-v1.0`.
