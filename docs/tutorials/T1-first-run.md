# T1. First Run in Five Minutes

| | |
| --- | --- |
| **Level** | Introductory |
| **Version** | 1.0 (2026-11-04) |
| **Time** | about 5 minutes, plus the build |
| **Timings measured on** | Apple M5 Pro (18 cores), 64 GB of RAM, macOS 26.6.2, Java 21.0.12 (Oracle JDK), Maven 3.9.16 |
| **Prerequisites** | Java 21, Maven 3.9+, git |

This tutorial builds Transept, runs it on the smallest dataset of the benchmark, and walks
through the files it writes. At the end you will know how to launch a run, where its results are,
how to read them, and how to repeat a run exactly.

## Step 1: get the code and the data

The tool and the benchmark data live in two repositories. Clone them side by side, and use the
published version of the dataset:

```bash
mkdir MSA && cd MSA
git clone https://github.com/jMetal/Transept-core.git
git clone https://github.com/TM-M2StructAlign/Dataset_Structural_GPCRs.git
git -C Dataset_Structural_GPCRs checkout v1.0.0
```

The commands of all tutorials are run from `MSA/Transept-core` and refer to the data as
`../Dataset_Structural_GPCRs/GPCRdb`.

## Step 2: build

```bash
cd Transept-core
mvn clean package -DskipTests
```

The result is `target/transept-core-0.9-jar-with-dependencies.jar`, a single file with every
dependency. `mvn clean verify` also runs the tests (a few minutes).

Maven downloads every dependency from Maven Central the first time, including
[jMetal](https://github.com/jMetal/jMetal) 7.11 and [Evolver](https://github.com/jMetal/Evolver)
2.6, so the first build takes longer.

## Step 3: the dataset

We use `classF`: the 11 human class F GPCRs (frizzled receptors and smoothened). A run reads four
kinds of files from `GPCRdb`:

| File | Content |
| --- | --- |
| `sequences/fasta/classF.fasta` | The 11 sequences |
| `sequences/tmregions/classF_predicted_topologies.3line` | Predicted topology: a line per residue type (`M` transmembrane, `I` inside, `O` outside, `S` signal peptide) |
| `distances/classF/<protein>_D.csv`, `<protein>_idx.csv` | Cα distance matrix of the AlphaFold2 model and the PDB residue numbers of its rows |
| `precomputed/classF/<tool>_classF.fasta` | Alignments produced by other tools, used as seeds of the search |

Transept does not align from scratch: it starts from the alignments of other tools (the
*seeds*) and searches for alignments that improve two objectives at the same time, a
topology-aware sum of pairs (**SOPwTP**) and a structural score (**lDDT**). Both are explained in
[T2](T2-problem-formulation.md).

## Step 4: run

```bash
java -jar target/transept-core-0.9-jar-with-dependencies.jar run \
  --data-root ../Dataset_Structural_GPCRs/GPCRdb \
  --dataset classF \
  --seeds clustalw,kalign,mafft,muscle5,tcoffee,tmaligner \
  --output-root results \
  --run-id first \
  --seed 1 \
  --max-evaluations 500 \
  --population-size 20
```

| Option | Meaning here |
| --- | --- |
| `--seeds` | Start from the alignments of six sequence aligners |
| `--output-root`, `--run-id` | Write to `results/ejecuciones/classF/first/` |
| `--seed` | Seed of the random generator |
| `--max-evaluations`, `--population-size` | A tiny run: about 500 alignments evaluated, a population of 20 (NSGA-II creates 100 offspring per generation, so it stops at 520) |

The algorithm is NSGA-II, the default. After some log lines from jMetal, the run ends with:

```text
Run finished: 520 evaluations, 20 solutions in results/ejecuciones/classF/first
```

It takes about 6 seconds. The run directory must not exist: to repeat a run, give it another
`--run-id`.

## Step 5: read the results

```text
results/ejecuciones/classF/first/
├── FUN.csv            objective values, one line per solution
├── MSASol0.fasta      alignment of line 0 of FUN.csv
├── ...
├── MSASol19.fasta
├── config.json        configuration of the run
├── metadata.json      provenance: commit, seeds, structures
├── progress.jsonl     progress of the search
├── runtime.csv        running time
└── seed.txt           random seed
```

**`FUN.csv`** has one line per solution, with SOPwTP and lDDT. Higher is better for both:

```text
60495.0,0.8001938417574468
58085.0,0.8125907400854824
60495.0,0.8001938417574468
60343.0,0.8043033089926177
58628.0,0.8103732752287465
...
```

The solutions trade one objective for the other: the first has the best SOPwTP and the worst
lDDT, the second the opposite. The file holds the final population, so some lines repeat (here
the first and the third): the search can produce copies of the
same alignment.

**`MSASol<i>.fasta`** is the alignment of line `i`, in aligned FASTA. Any alignment viewer (Jalview,
AliView, ...) opens it:

```text
>fzd2_human
MR----PRSALP-----------------------------------------RLLL---
----PLLLL-----PAAGPAQ--------------------FHGEKGISIPD--------
...
```

**`config.json`** records every option with its default resolved, and **`metadata.json`** what is
needed to trust and repeat the run: the git commit of the build, the SHA-256 of every seed file,
and how each AlphaFold model was mapped to its sequence. In `classF` every protein has the status
`NUMBERED` (the residue numbers of the model match the sequence) and `proteinsWithoutStructure` is
empty. Look at these two fields whenever you use a new dataset: a protein without structure does
not contribute to the lDDT.

**`progress.jsonl`** has one line every 100 evaluations with the best value of each objective so far.
It can be followed while a long run is active (`tail -f`):

```json
{"evaluations":120,"elapsedMs":1692,"best":{"SOPwTP":60495.0,"lDDT":0.8075674397684053}}
```

The [output reference](../reference/outputs.md) describes every field.

## Step 6: repeat a run exactly

The same options and the same seed give the same result, bit for bit:

```bash
java -jar target/transept-core-0.9-jar-with-dependencies.jar run \
  --data-root ../Dataset_Structural_GPCRs/GPCRdb --dataset classF \
  --seeds clustalw,kalign,mafft,muscle5,tcoffee,tmaligner \
  --output-root results --run-id first-again --seed 1 \
  --max-evaluations 500 --population-size 20

cmp results/ejecuciones/classF/first/FUN.csv results/ejecuciones/classF/first-again/FUN.csv \
  && echo identical
```

This holds for the same version of the tool, whatever the number of `--threads`. A different
`--seed` gives a different, equally valid result: the search is stochastic, and a study runs
several seeds ([T3](T3-several-seeds.md)).

## Step 7: a real run

A useful run needs more evaluations and a larger population. With the defaults of the command
line (10 000 evaluations, populations of 100) and four evaluation threads it takes about
30 seconds on `classF`:

```bash
java -jar target/transept-core-0.9-jar-with-dependencies.jar run \
  --data-root ../Dataset_Structural_GPCRs/GPCRdb --dataset classF \
  --seeds clustalw,kalign,mafft,muscle5,tcoffee,tmaligner \
  --output-root results --run-id seed_000 --seed 0 \
  --max-evaluations 10000 --population-size 100 --threads 4
```

Larger datasets take much longer: the cost grows with the square of the number of proteins and
with their length, so `classA_002` (76 proteins) is about two orders of magnitude slower.

## What next

- [T2](T2-problem-formulation.md) explains what the algorithm optimises and how it changes the
  alignments.
- [T4](T4-choosing-and-evaluating.md) shows how to choose alignments from the result and compare
  them with the GPCRdb reference alignments.
- The [command-line reference](../reference/cli.md) lists every option.
