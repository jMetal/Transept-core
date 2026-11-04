# T5. Comparing Algorithms and Operators

| | |
| --- | --- |
| **Level** | Intermediate |
| **Version** | 1.0 (2026-11-04) |
| **Time** | about 20 minutes, plus about 10 minutes of runs |
| **Timings measured on** | Apple M5 Pro (18 cores), 64 GB of RAM, macOS 26.6.2, Java 21.0.12 (Oracle JDK), Python 3.9 with NumPy and Matplotlib 3.9 |
| **Prerequisites** | [T3. A study with several seeds](T3-several-seeds.md) |

TM-M2StructAlign can run several multiobjective algorithms and operators on the same problem. This
tutorial compares five configurations on `classF` with a reduced budget: how to run them so that
the comparison is fair, how to measure them, and how to tell a real difference from noise.

## Step 1: what can be compared

| Option | Values |
| --- | --- |
| `--algorithm` | `NSGA_II` (default), `NSGA_III`, `AGE_MOEA`, `AGE_MOEA_II`, `SMS_EMOA`, `MOEA_D`, `PAES`, and the random searches `RANDOM_SEARCH` and `RANDOM_SEEDED` |
| `--mutation` | `GAP_SHIFT` (default), `TM_SAFE_GAP_SHIFT`, `MERGE_GAP_RUNS`, `RANDOM_CHOICE` |
| `--crossover` | `SINGLE_POINT` (default), `TWO_POINT` |
| `--initialization` | `CROSSOVER` (default), `CROSSOVER_DISTINCT`, `PERTURBED` |

All algorithms start from the same seed alignments, use the same operators and stop after the same
number of evaluations, so that only the option under study changes. Some need care:

- `NSGA_III` needs a population size that is a multiple of 4.
- `SMS_EMOA` and `MOEA_D` are steady state: they evaluate one solution at a time, so `--threads`
  does not speed them up.
- `PAES` is a (1+1) strategy: it starts from the first seed and only mutates (mutation probability
  1.0 by default).
- The random searches are **sanity checks**. `RANDOM_SEEDED` applies 10 random mutations to copies of
  the seeds: an algorithm that does not beat it is not using its operators well. `RANDOM_SEARCH`
  places gaps at random and ignores the seeds: a much weaker baseline.

## Step 2: run the configurations

Name every run `<configuration>_seed_<seed>`, so that the analysis can group and pair them, and
use the same seeds for every configuration:

```bash
for seed in 1 2 3 4 5 6; do
  for config in "nsgaii:--algorithm NSGA_II" \
                "nsgaii-tmsafe:--algorithm NSGA_II --mutation TM_SAFE_GAP_SHIFT" \
                "smsemoa:--algorithm SMS_EMOA" \
                "moead:--algorithm MOEA_D" \
                "randomseeded:--algorithm RANDOM_SEEDED"; do
    name=${config%%:*}; options=${config#*:}
    java -jar target/TM-MSAStructAligner-1.0-jar-with-dependencies.jar run \
      --data-root ../Dataset_Structural_GPCRs/GPCRdb --dataset classF \
      --seeds clustalw,kalign,mafft,muscle5,tcoffee,tmaligner \
      --output-root comparison --run-id ${name}_seed_$seed --seed $seed \
      --max-evaluations 2500 --population-size 100 --threads 4 $options
  done
done
```

The 30 runs take about 10 minutes. A run of NSGA-II takes about 8 seconds; SMS-EMOA and MOEA/D
about 32 seconds, because they cannot evaluate in parallel.

## Step 3: measure every run

Analyse all the runs together, so that the objectives are normalised over the same range for every
configuration:

```bash
python3 scripts/pareto_analysis.py comparison classF --output comparison/analysis_classF
```

```text
Analyzed 2425 solutions from 30 runs
Combined non-dominated objective vectors: 51
```

`indicators_by_run.csv` has the hypervolume and the IGD+ of every run ([T3](T3-several-seeds.md),
Step 3).

> **Keep the weak baselines out of the joint analysis.** The normalisation uses the range of all
> the solutions analysed. `RANDOM_SEARCH` produces alignments far worse than the rest, which would
> stretch that range and squeeze the differences between the other configurations. Analyse it
> separately if you need it.

## Step 4: compare

`scripts/compare_runs.py` groups the runs by configuration, pairs them by seed, and compares every
configuration with a baseline:

```bash
python3 scripts/compare_runs.py comparison/analysis_classF/indicators_by_run.csv --baseline nsgaii
```

```text
hypervolume (higher is better), baseline nsgaii
configuration        runs   median      IQR    wins       p
nsgaii                  6    0.841    0.130
moead                   6    0.741    0.052     1/6   0.094
nsgaii-tmsafe           6    0.809    0.134     3/6   1.000
randomseeded            6    0.303    0.116     0/6   0.031
smsemoa                 6    0.811    0.115     2/6   0.688
```

`--indicator igd_plus` compares by IGD+ instead:

```text
igd_plus (lower is better), baseline nsgaii
configuration        runs   median      IQR    wins       p
nsgaii                  6    0.026    0.028
moead                   6    0.056    0.024     0/6   0.031
nsgaii-tmsafe           6    0.026    0.021     3/6   1.000
randomseeded            6    0.443    0.142     0/6   0.031
smsemoa                 6    0.026    0.022     3/6   1.000
```

| Column | Meaning |
| --- | --- |
| `median`, `IQR` | Median and interquartile range of the indicator over the seeds |
| `wins` | Seeds on which the configuration is better than the baseline with the same seed |
| `p` | Two-sided p-value of the exact Wilcoxon signed-rank test on the paired runs |

## Step 5: interpret

- **The search works.** Every algorithm is far better than `RANDOM_SEEDED` (hypervolume 0.30
  against 0.74–0.84): the operators and the selection add much more than random perturbations of
  the seeds. With six seeds, 0.031 is the smallest possible p-value, reached when one
  configuration wins on every seed.
- **MOEA/D is worse than NSGA-II** on this problem: it loses on every seed in IGD+ (p = 0.031) and
  contributes 1 of the 51 solutions of the combined front.
- **NSGA-II, NSGA-II with the TM-safe mutation and SMS-EMOA are indistinguishable** with these
  runs: the differences in median are smaller than the spread, and the p-values are large.
  SMS-EMOA contributes the most solutions to the combined front (69 against 46 of NSGA-II), but
  costs four times as much time.
- **"Not significant" does not mean "equal".** Six seeds detect only large differences. To claim
  that two configurations are equivalent, or to detect small differences, use more seeds (10 to
  30), more than one dataset, and a budget close to the one of the real study.

## Step 6: good practice

- **Change one thing at a time**, with everything else at its default, and keep the same seeds for
  every configuration: the test pairs runs by seed.
- **Report the spread**, not only the medians, and the number of seeds.
- **Use several datasets.** The result on one dataset does not always carry over: while the tool
  was developed, some operators made a significant difference on `classF` and none on
  `classT2_19`.
- **Look beyond the indicators.** The hypervolume measures the objectives, which are not the whole
  story ([T4](T4-choosing-and-evaluating.md)): compare also the chosen alignments with the reference
  alignments when they exist.
- **Account for the cost.** Steady-state algorithms evaluate one solution at a time and do not use
  `--threads`. At equal evaluations they cost more time.
