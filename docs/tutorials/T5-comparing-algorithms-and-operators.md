# T5. Comparing Algorithms and Operators

| | |
| --- | --- |
| **Level** | Intermediate |
| **Version** | 1.0 (2026-11-04) |
| **Time** | about 20 minutes, plus about 10 minutes of runs |
| **Timings measured on** | Apple M5 Pro (18 cores), 64 GB of RAM, macOS 26.6.2, Java 21.0.12 (Oracle JDK), Python 3.9 with NumPy and Matplotlib 3.9 |
| **Prerequisites** | [T3. A study with several seeds](T3-several-seeds.md) |

Transept can run several multiobjective algorithms and operators on the same problem. This
tutorial compares five configurations on `classF` with a reduced budget: how to run them so that
the comparison is fair, how to measure them, and how to tell a real difference from noise.

## Step 1: what can be compared

| Option | Values |
| --- | --- |
| `--algorithm` | `NSGA_II` (default), `NSGA_III`, `AGE_MOEA`, `AGE_MOEA_II`, `SMS_EMOA`, `MOEA_D`, `PAES`, and the random searches `RANDOM_SEARCH` and `RANDOM_SEEDED` |
| `--mutation` | `TM_SAFE_GAP_SHIFT` (default), `MERGE_GAP_RUNS`, `RANDOM_CHOICE` |
| `--objectives` | `SOPWTP_LDDT` (default), `LDDT` (with `--algorithm GENETIC_ALGORITHM`), `SOPWTP_LDDT_COVERAGE` |
| `--configuration` | NSGA-II only: any configuration of its parameter space, e.g. `--selection random` ([configuration](../reference/configuration.md)) |

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
                "nsgaii-randomchoice:--algorithm NSGA_II --mutation RANDOM_CHOICE" \
                "smsemoa:--algorithm SMS_EMOA" \
                "moead:--algorithm MOEA_D" \
                "randomseeded:--algorithm RANDOM_SEEDED"; do
    name=${config%%:*}; options=${config#*:}
    java -jar target/transept-core-0.1.0-jar-with-dependencies.jar run \
      --data-root ../Dataset_Structural_GPCRs/GPCRdb --dataset classF \
      --seeds clustalw,kalign,mafft,muscle5,tcoffee,tmaligner \
      --output-root comparison --run-id ${name}_seed_$seed --seed $seed \
      --max-evaluations 2500 --population-size 100 --threads 4 $options
  done
done
```

The configurations are the default (NSGA-II with the `TM_SAFE_GAP_SHIFT` mutation), NSGA-II with
the `RANDOM_CHOICE` mutation, two other algorithms and the strong sanity check. The 30 runs take
about 9 minutes. A run of NSGA-II takes about 9 seconds; SMS-EMOA and MOEA/D about 31 seconds,
because they cannot evaluate in parallel.

## Step 3: measure every run

Analyse all the runs together, so that the objectives are normalised over the same range for every
configuration:

```bash
python3 scripts/pareto_analysis.py comparison classF --output comparison/analysis_classF
```

```text
Analyzed 2410 solutions from 30 runs
Combined non-dominated objective vectors: 88
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
nsgaii                  6    0.925    0.035
moead                   6    0.898    0.030     0/6   0.031
nsgaii-randomchoice     6    0.939    0.023     5/6   0.438
randomseeded            6    0.032    0.051     0/6   0.031
smsemoa                 6    0.935    0.013     4/6   1.000
```

`--indicator igd_plus` compares by IGD+ instead:

```text
igd_plus (lower is better), baseline nsgaii
configuration        runs   median      IQR    wins       p
nsgaii                  6    0.019    0.013
moead                   6    0.033    0.009     0/6   0.031
nsgaii-randomchoice     6    0.011    0.008     5/6   0.062
randomseeded            6    1.032    0.169     0/6   0.031
smsemoa                 6    0.018    0.007     4/6   1.000
```

| Column | Meaning |
| --- | --- |
| `median`, `IQR` | Median and interquartile range of the indicator over the seeds |
| `wins` | Seeds on which the configuration is better than the baseline with the same seed |
| `p` | Two-sided p-value of the exact Wilcoxon signed-rank test on the paired runs |

## Step 5: interpret

- **The search works.** Every algorithm is far better than `RANDOM_SEEDED` (hypervolume 0.03
  against 0.90–0.94): the operators and the selection add much more than random perturbations of
  the seeds. With six seeds, 0.031 is the smallest possible p-value, reached when one
  configuration wins on every seed.
- **MOEA/D is worse than NSGA-II** on this problem: it loses on every seed in both indicators
  (p = 0.031) and contributes none of the 88 solutions of the combined front.
- **`RANDOM_CHOICE` looks promising but is not proven.** It wins on 5 of 6 seeds and contributes
  most of the combined front (111 solutions, against 20 of NSGA-II with the default mutation), but
  the p-values (0.44 and 0.06) do not support a difference. It is a candidate for a larger study,
  not a reason to change the default.
- **NSGA-II and SMS-EMOA are indistinguishable** with these runs: the differences in median are
  smaller than the spread, and the p-values are large. SMS-EMOA costs three times as much time.
- **"Not significant" does not mean "equal".** Six seeds detect only large differences. To claim
  that two configurations are equivalent, or to detect small differences, use more seeds (10 to
  30), more than one dataset, and a budget close to the one of the real study.

Formulations with a different number of objectives cannot be compared by hypervolume: each one
measures a different space. Compare them instead by the agreement with the reference alignments
of an alignment chosen with a rule fixed in advance (for example, the best lDDT), and by the
spread of that agreement over the returned alignments ([T4](T4-choosing-and-evaluating.md)).
`scripts/formulation_study.py` does this for every run of a dataset (rule R1: best lDDT; R2: the
knee) and writes a file that `compare_runs.py` reads:

```bash
python3 scripts/formulation_study.py comparison classF \
  --data-root ../Dataset_Structural_GPCRs/GPCRdb --output comparison/formulations_classF.csv
python3 scripts/compare_runs.py comparison/formulations_classF.csv --baseline nsgaii \
  --indicator r1_pair_f1
```

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
