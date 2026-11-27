# T8. Repairing the Gaps Inside Helices

| | |
| --- | --- |
| **Level** | Intermediate |
| **Version** | 1.0 (2026-11-19) |
| **Time** | about 15 minutes |
| **Prerequisites** | [T1. First run in five minutes](T1-first-run.md) |
| **Checked with** | The 11 aligners of the GPCR benchmark on 15 datasets (Step 5) |

Most aligners do not know where the transmembrane helices are, and place gaps inside them: a gap
inside a helix means that a stretch of a helix has no counterpart in the other proteins, which is
rarely what the structures say. Transept's operators never open such a gap, but the seed alignments
bring their own. The **repair** moves them out, without any search: before every run (unless
`--no-repair-seeds` is given), or on any alignment with the `repair` command.

This tutorial explains what the repair does, how to use it, and how much it changes the accuracy of
the aligners of the benchmark.

## Step 1: what counts as a gap inside a helix

In every row, a gap run *splits a helix* when it lies between two consecutive residues of the same
transmembrane helix (`M` in the topology). A gap run before the first residue of a helix or after
its last residue does not split it. Gaps in the other rows do not matter here: the repair looks at
one protein at a time.

## Step 2: what the repair does

For every row, the repair finds each gap run that splits a helix and moves it, whole, to the nearest
position that does not split one, counting residues; on a tie it goes towards the N-terminus. It
repeats this until no run splits a helix. Three things never change:

- the residues of every row keep their order;
- the length of the alignment does not change: moving a run of `k` gaps across `d` residues shifts
  those `d` residues `k` columns;
- rows without gaps inside helices are left as they were.

A real case, from the MAFFT alignment of `classF` in the benchmark. MAFFT puts two gaps right after
the first residue of a helix of `fzd4_human`; the nearest boundary is one residue back, before that
residue:

```text
                 before                    after
residues    QEGLKNT--GCAIIF           QEGLKN--TGCAIIF
topology    OOOOOOM  MMMMMM           OOOOOO  MMMMMMM
```

The helix is now whole: `T` moved two columns to the right, next to the rest of its helix.

## Step 3: repairing alignments with the `repair` command

```bash
java -jar target/transept-core-1.0-SNAPSHOT-jar-with-dependencies.jar repair \
  --topology ../Dataset_Structural_GPCRs/GPCRdb/sequences/tmregions/classF_predicted_topologies.3line \
  --output-dir repaired \
  ../Dataset_Structural_GPCRs/GPCRdb/precomputed/classF/*_classF.fasta
```

It reads every alignment, checks that its rows are the proteins of the topology (with the same
residues), and writes the repaired alignment to `repaired/` with the same file name, reporting how
many gap runs it moved:

```text
mafft_classF.fasta: 4 gap runs moved -> repaired/mafft_classF.fasta
```

An alignment of other proteins, or with other residues, is an error; the command never overwrites
its inputs. To score the repaired alignments against a reference, use `score` with
`--msa-root repaired --glob "*.fasta"` ([command line](../reference/cli.md#score)).

## Step 4: repairing the seeds of a run

`run` applies the same repair to every seed alignment before the search. With `--no-repair-seeds` the
precomputed alignments keep their gaps, and so may the alignments of the front; the internal seeds
are always repaired, because their aligner does not know the topology. `metadata.json` records how many gap runs it moved in each seed (`repairedGapRuns`).
Since the operators never open a gap inside a helix, no alignment of the final front has one.

The graphical interface of Transept repairs the seeds by default (**Repair the seeds** in **New
run**).

## Step 5: how much the repair changes an alignment

The 11 aligners of the GPCR benchmark, repaired without search and scored against the GPCRdb
references (median over the datasets where they can be scored):

| Aligner | Datasets | Gap runs moved | Pair F1 before | Pair F1 change |
| --- | --- | --- | --- | --- |
| ClustalW | 15 | 27 | 0.737 | +0.002 |
| MAFFT | 15 | 51 | 0.691 | +0.004 |
| Kalign | 15 | 77 | 0.672 | +0.001 |
| MUSCLE5 | 15 | 79 | 0.730 | −0.003 |
| T-Coffee | 15 | 58 | 0.744 | −0.005 |
| TM-Aligner | 15 | 45 | 0.691 | −0.001 |
| FoldMason | 11 | 17 | 0.724 | +0.001 |
| US-align | 11 | 23 | 0.784 | −0.003 |
| MUSTANG | 11 | 25 | 0.774 | −0.001 |
| mTM-align | 11 | 65 | 0.798 | −0.022 |
| Caretta | 11 | 15 | 0.175 | +0.000 |

- After the repair, no alignment has a gap inside a helix.
- The accuracy hardly changes: within ±0.005 for 10 of the 11 aligners. The repair is a local move,
  not a realignment, so it neither fixes nor spoils much.
- mTM-align is the exception (−0.022): it leaves many gaps inside helices, and moving them probably
  shifts residues that it had aligned as the reference does.
- The TM-gap rate ([metrics](../reference/metrics.md#tm-gap-rate)) can go up after the repair: it
  counts transmembrane residues facing a gap in *another* protein, and a repaired helix can end up
  facing the gaps of other rows. Gaps inside helices and the TM-gap rate measure different things.

In short: on these datasets the repair removes the gaps inside helices at almost no cost in accuracy.
It is not what makes a Transept alignment better than its seeds; the search is. In more divergent
families the cost can be higher, because some gaps inside annotated helices are also in the reference
alignments: compare both options if you are not sure.

## Summary

| To | Use |
| --- | --- |
| Remove the gaps inside helices of any alignment | `repair --topology ... --output-dir ... ALIGNMENT...` |
| Start a run from repaired seeds | `run` (default) |
| Keep the gaps of the precomputed seeds | `run --no-repair-seeds` |
| See how many gap runs were moved in each seed | `repairedGapRuns` in `metadata.json` |
