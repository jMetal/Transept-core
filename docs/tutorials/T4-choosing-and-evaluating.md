# T4. Choosing an Alignment and Evaluating It Against GPCRdb

| | |
| --- | --- |
| **Level** | Intermediate |
| **Version** | 1.0 (2026-11-04) |
| **Time** | about 20 minutes, plus one run of 30 seconds |
| **Timings measured on** | Apple M5 Pro (18 cores), 64 GB of RAM, macOS 26.6.2, Java 21.0.12 (Oracle JDK), Python 3.9 |
| **Prerequisites** | [T1. First run in five minutes](T1-first-run.md), [T2. Problem formulation](T2-problem-formulation.md), Python 3 |

A run returns many alignments, each a different compromise between the two objectives. This
tutorial shows how to inspect them, how to pick one, and how to check the choice against the
reference alignments of GPCRdb, which the search never sees.

## Step 1: a run

```bash
java -jar target/transept-core-0.1.0-jar-with-dependencies.jar run \
  --data-root ../Dataset_Structural_GPCRs/GPCRdb --dataset classF \
  --seeds clustalw,kalign,mafft,muscle5,tcoffee,tmaligner \
  --output-root results --run-id seed_042 --seed 42 \
  --max-evaluations 10000 --population-size 100 --threads 4
```

This is the default configuration (NSGA-II, `TM_SAFE_GAP_SHIFT` mutation), with seed 42. It takes
about 30 seconds.

## Step 2: the non-dominated alignments

`FUN.csv` holds the final population, which can contain dominated and repeated solutions
([outputs](../reference/outputs.md)). `scripts/summarize_front.py` keeps the distinct
non-dominated ones and sorts them from the best SOPwTP to the best lDDT:

```bash
python3 scripts/summarize_front.py results/ejecuciones/classF/seed_042
```

```text
100 solutions, 88 distinct non-dominated
file            SOPwTP    lDDT    ...  note
MSASol0.fasta     61332  0.8002        best SOPwTP
MSASol78.fasta    61267  0.8004
...
MSASol47.fasta    59498  0.8114        knee
...
MSASol92.fasta    56196  0.8207
MSASol1.fasta     56171  0.8209        best lDDT
```

Here every solution of the population is non-dominated, but 12 of them repeat another one. Three
alignments stand out:

- the **extremes**: the best SOPwTP (`MSASol0`) and the best lDDT (`MSASol1`);
- the **knee** (`MSASol47`): after scaling both objectives to [0, 1] over the front, the solution
  closest to the ideal point (1, 1). It is a common default when no objective is preferred.

The front is smooth: along it, every gain of lDDT costs some SOPwTP, without jumps. A jump (a
large loss of one objective for a tiny gain of the other) would deserve a closer look before
trusting the alignments beyond it.

## Step 3: compare with the reference alignments

GPCRdb provides structure-based reference alignments of the benchmark datasets. `EvaluationMain`
compares every alignment of a run with the reference ([metrics](../reference/metrics.md)):

```bash
java -cp target/transept-core-0.1.0-jar-with-dependencies.jar \
  org.transept.evaluation.EvaluationMain \
  --reference ../Dataset_Structural_GPCRs/GPCRdb/reference_alignments/classF.fasta \
  --topology ../Dataset_Structural_GPCRs/GPCRdb/sequences/tmregions/classF_predicted_topologies.3line \
  --msa-root results/ejecuciones/classF/seed_042 \
  --output results/evaluation_seed_042.csv
```

| Metric | Meaning | Better |
| --- | --- | --- |
| Pair F1 | Agreement of the aligned residue pairs with the reference | Higher |
| TM-pair F1 | The same, restricted to pairs of transmembrane residues | Higher |
| TM-gap rate | How often a transmembrane residue faces a gap in another protein | Lower |

Give the evaluation file to the summary script to see both views side by side:

```bash
python3 scripts/summarize_front.py results/ejecuciones/classF/seed_042 \
  results/evaluation_seed_042.csv
```

```text
file            SOPwTP    lDDT    pair-F1  TM-pair-F1  TM-gap-rate  note
MSASol0.fasta     61332  0.8002    0.5380      0.9995       0.0061  best SOPwTP
...
MSASol47.fasta    59498  0.8114    0.5555      0.9998       0.0052  knee
...
MSASol92.fasta    56196  0.8207    0.5637      0.9998       0.0044
MSASol1.fasta     56171  0.8209    0.5639      0.9998       0.0045  best lDDT
```

## Step 4: interpret

- **The pair F1 grows with the lDDT.** Along the front, from the best SOPwTP to the best lDDT, the
  agreement with the reference goes from 0.538 to 0.564, almost monotonically. Here the
  best-lDDT extreme is the best alignment of the run, and the knee lies in the middle (0.5555).
- **No alignment breaks a helix.** The TM-gap rate stays between 0.0044 and 0.0061 and comes from
  the seeds: neither the crossover nor the mutation ever puts a gap inside a transmembrane helix
  ([T2](T2-problem-formulation.md)).
- **Transmembrane regions are easy in classF.** Every alignment has a TM-pair F1 of at least
  0.9995: the differences are in the loops and the termini.
- **A local structural score is not enough.** The lDDT is superposition-free and local: it hardly
  notices a helix out of register. That it agrees with the reference here does not make it a
  substitute for one; check the chosen alignments whenever a reference exists.

## Step 5: compare with the seeds

The search starts from the seeds, so the first question is whether it improves them. Evaluate the
seed alignments with the same command. Use `--glob` to select them, and always write the report
outside the dataset directory:

```bash
java -cp target/transept-core-0.1.0-jar-with-dependencies.jar \
  org.transept.evaluation.EvaluationMain \
  --reference ../Dataset_Structural_GPCRs/GPCRdb/reference_alignments/classF.fasta \
  --topology ../Dataset_Structural_GPCRs/GPCRdb/sequences/tmregions/classF_predicted_topologies.3line \
  --msa-root ../Dataset_Structural_GPCRs/GPCRdb/precomputed/classF \
  --glob '{clustalw,kalign,mafft,muscle5,tcoffee,tmaligner}_classF.fasta' \
  --output results/evaluation_seeds_classF.csv
```

| Seed | Pair F1 | TM-pair F1 | TM-gap rate |
| --- | --- | --- | --- |
| ClustalW | 0.5598 | 0.9998 | 0.0043 |
| MUSCLE5 | 0.5515 | 0.9998 | 0.0022 |
| T-Coffee | 0.5511 | 0.9998 | 0.0018 |
| TM-Aligner | 0.5463 | 0.9983 | 0.0036 |
| MAFFT | 0.5442 | 0.9998 | 0.0052 |
| Kalign | 0.5413 | 0.9995 | 0.0032 |

The knee (pair F1 0.5555) is better than five of the six seeds; the best-lDDT extreme (0.5639) is
better than all of them, ClustalW (0.5598) included. The gains are small: on this dataset the
search mostly recombines and polishes what the seeds already have.

## Step 6: draw conclusions carefully

- **One run is one sample.** The search is stochastic: repeat it with several seeds and compare
  distributions ([T3](T3-several-seeds.md), [T5](T5-comparing-algorithms-and-operators.md)).
- **Choose the rule before looking at the reference.** Picking the alignment with the best pair F1
  uses the reference to choose and overstates what the method achieves. Report the knee or the
  best-lDDT extreme, or the rule you decided in advance, and use the reference only to evaluate.
- **The reference exists only for the benchmark.** On your own proteins
  ([T6](T6-your-own-dataset.md)), you choose with the objectives alone. In the benchmark the
  best-lDDT extreme agreed best with the reference in most of the runs made while developing the
  tool, but the knee is a safer default when the two objectives disagree strongly.
