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
java -jar target/TM-MSAStructAligner-1.0-jar-with-dependencies.jar run \
  --data-root ../Dataset_Structural_GPCRs/GPCRdb --dataset classF \
  --seeds clustalw,kalign,mafft,muscle5,tcoffee,tmaligner \
  --output-root results --run-id seed_042 --seed 42 \
  --max-evaluations 10000 --population-size 100 --threads 4
```

This is the default configuration (NSGA-II, `GAP_SHIFT` mutation), with seed 42. It takes about
30 seconds.

## Step 2: the non-dominated alignments

`FUN.csv` holds the final population, which can contain dominated and repeated solutions
([outputs](../reference/outputs.md)). `scripts/summarize_front.py` keeps the distinct
non-dominated ones and sorts them from the best SOPwTP to the best lDDT:

```bash
python3 scripts/summarize_front.py results/ejecuciones/classF/seed_042
```

```text
100 solutions, 81 distinct non-dominated
file            SOPwTP    lDDT    ...  note
MSASol0.fasta     61113  0.5651        best SOPwTP
MSASol97.fasta    61084  0.5658
...
MSASol6.fasta     58204  0.6062        knee
MSASol3.fasta     57106  0.6072
MSASol4.fasta     50107  0.6074
...
MSASol1.fasta     48650  0.6141        best lDDT
```

Here every solution of the population is non-dominated, but 19 of them repeat another one. Three
alignments stand out:

- the **extremes**: the best SOPwTP (`MSASol0`) and the best lDDT (`MSASol1`);
- the **knee** (`MSASol6`): after scaling both objectives to [0, 1] over the front, the solution
  closest to the ideal point (1, 1). It is a common default when no objective is preferred.

Look at the values between `MSASol3` and `MSASol4`: the front has a **jump**. The last seven
alignments gain less than 0.007 of lDDT and lose about 7 000 of SOPwTP. A jump like this deserves
a closer look before trusting the alignments beyond it.

## Step 3: compare with the reference alignments

GPCRdb provides structure-based reference alignments of the benchmark datasets. `EvaluationMain`
compares every alignment of a run with the reference ([metrics](../reference/metrics.md)):

```bash
java -cp target/TM-MSAStructAligner-1.0-jar-with-dependencies.jar \
  org.tmm2structalign.evaluation.EvaluationMain \
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
MSASol0.fasta     61113  0.5651    0.5448      0.9998       0.0052  best SOPwTP
...
MSASol6.fasta     58204  0.6062    0.5566      0.9998       0.0025  knee
MSASol3.fasta     57106  0.6072    0.5557      0.9998       0.0043
MSASol4.fasta     50107  0.6074    0.5144      0.9354       0.0624
...
MSASol1.fasta     48650  0.6141    0.5143      0.9354       0.0648  best lDDT
```

## Step 4: interpret

- **The knee is the best choice here.** It shares the highest pair F1 of the front (0.5566), and half
  the TM-gap rate of the best-SOPwTP extreme.
- **The best-lDDT extreme is the worst alignment.** Beyond the jump, the alignments put gaps inside
  transmembrane helices: the TM-gap rate grows from 0.0025 to 0.062 and the TM-pair F1 falls from
  0.9998 to 0.935. The default mutation only *discourages* such gaps, and the lDDT rewards them a
  little because it favours compact alignments and hardly notices a helix out of register
  ([T2](T2-problem-formulation.md)). Do not choose by the second objective alone.
- **Transmembrane regions are easy in classF.** Every alignment except those beyond the jump has a
  TM-pair F1 of 0.9998: the differences are in the loops and the termini.

## Step 5: compare with the seeds

The search starts from the seeds, so the first question is whether it improves them. Evaluate the
seed alignments with the same command. Use `--glob` to select them, and always write the report
outside the dataset directory:

```bash
java -cp target/TM-MSAStructAligner-1.0-jar-with-dependencies.jar \
  org.tmm2structalign.evaluation.EvaluationMain \
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

The knee (pair F1 0.5566) is better than five of the six seeds, but not than ClustalW (0.5598).
On this dataset and with this budget, the search mostly recombines what the seeds already have.

## Step 6: draw conclusions carefully

- **One run is one sample.** The search is stochastic: repeat it with several seeds and compare
  distributions ([T3](T3-reproducing-the-paper.md), [T5](T5-comparing-algorithms-and-operators.md)).
- **Choose the rule before looking at the reference.** Picking the alignment with the best pair F1
  uses the reference to choose and overstates what the method achieves. Report the knee, or the
  rule you decided in advance, and use the reference only to evaluate.
- **The reference exists only for the benchmark.** On your own proteins
  ([T6](T6-your-own-dataset.md)), you choose with the objectives alone: prefer the knee, and look at
  the TM-gap behaviour of the extremes with the topology in hand.
