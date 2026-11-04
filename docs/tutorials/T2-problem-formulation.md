# T2. Problem Formulation, Encoding and Variation Operators

| | |
| --- | --- |
| **Level** | Introductory |
| **Version** | 1.0 (2026-11-04) |
| **Time** | about 45 minutes |
| **Prerequisites** | [T1. First run in five minutes](T1-first-run.md) |
| **Companion** | `org.tmm2structalign.example.FormulationTour` |

This tutorial describes what TM-M2StructAlign optimises and how it does it: the inputs, the two
objectives, the constraints, how an alignment is encoded, how the initial population is built,
and what every crossover and mutation does. A small example with three proteins runs through it.
Its output is printed by a class of the project, so you can run it and change it:

```bash
java -cp target/TM-MSAStructAligner-1.0-jar-with-dependencies.jar \
  org.tmm2structalign.example.FormulationTour
```

The output shown below is the one of this version; a test checks that the class still prints it.

## Step 1: the problem

**Input.** A set of N protein sequences and, for every protein:

- its predicted **topology**, one letter per residue: `M` (transmembrane), `I` (inside), `O`
  (outside) or `S` (signal peptide). Only the distinction between `M` and the rest matters;
- the **Cα distance matrix** of its AlphaFold2 model, with the residue numbers of the model so that
  the matrix can be mapped onto the sequence;
- a few **seed alignments** of the N sequences, produced by other tools.

**Output.** A set of multiple alignments of the N sequences that are good compromises between the
two objectives below: none of them is better than another in both.

The three proteins of the example have one transmembrane helix each. Their structures are
synthetic: residues that the seed alignment places in the same column sit at nearly the same
position of an ideal helix, so the seed alignment is the structurally correct one.

```text
== 1. Proteins ==
p1  PQTLLVAVLGNKE
    OOOMMMMMMMIII
p2  PTLLIAVLGKE
    OOMMMMMMMII
p3  PQSTLVAILGNE
    OOOOMMMMMMII
```

## Step 2: encoding

An alignment is stored as one **gap layout** per protein: the aligned length and the runs of gap
columns of its row (class `alignment.GapIntervals`). The residues are never stored: the aligned row
is decoded by inserting the gaps into the ungapped sequence.

```text
== 2. Alignment and encoding ==
p1  PQTLLVAVLGNKE  gap runs []
p2  P-TLLIAVLG-KE  gap runs [1, 10]
p3  PQSTLVAILGN-E  gap runs [11]
length 13
transmembrane columns (M), gaps included:
p1  ...MMMMMMM...
p2  ...MMMMMMM...
p3  ....MMMMMM...
```

The row of `p2` is encoded as length 13 and runs `[1, 1]` and `[10, 10]` (printed as `1` and `10`).
The encoding has three properties that the operators rely on:

- **Immutable and canonical.** Runs are sorted, inclusive and separated by at least one residue, so
  two rows with the same gaps are equal. Operators never modify a layout: they build a new one.
- **Valid by construction.** A layout always decodes to the right sequence; the only condition
  between rows is that all of them have the same length.
- **Small.** It takes memory proportional to the number of gap runs, not to the length of the
  alignment.

A solution (class `alignment.MsaSolution`) is the list of the N layouts in a fixed protein order,
plus the two objective values.

**Gaps also have a topology.** The last block of the output marks the transmembrane columns of
every row. A residue is transmembrane if its topology letter is `M`. A **gap is transmembrane when
the residues on both sides of it are transmembrane**: the gap of `p2` at column 1 lies between
`P` (outside) and `T` (outside), so it is not transmembrane, and a gap between two helix residues
would be. This rule decides the gap penalties of the first objective and what "breaking a helix"
means for the operators.

## Step 3: objective 1, SOPwTP

The **topology-aware sum of pairs** adds a score over every pair of sequences and every column of
the alignment:

| Column of the pair | Score |
| --- | --- |
| Two residues, both transmembrane | PHAT matrix (built from transmembrane segments) |
| Two residues, otherwise | BLOSUM62 |
| A gap and a residue | Affine gap penalty, subtracted: transmembrane (open 8, extend 3) if the gap or the residue is transmembrane, otherwise general (open 3, extend 1) |
| Two gaps | Skipped; it does not close an open gap |

A gap *extends* when the previous scored column of the pair had a gap in the same sequence, and
*opens* otherwise. The total is the sum over all pairs. Higher is better.

A worked example with two sequences, `a = KLIVG` (topology `OMMMO`) and `b = KLAIVSG` (topology
`OMMMMOO`):

```text
a  K L - I V - G
b  K L A I V S G
   5 4 -8 5 4 -3 6   = 13
```

| Column | Pair | Rule | Score |
| --- | --- | --- | --- |
| 0 | K/K | both outside: BLOSUM62 | +5 |
| 1 | L/L | both transmembrane: PHAT | +4 |
| 2 | gap/A | the gap lies between L and I, both transmembrane, so it is transmembrane: open 8 | −8 |
| 3 | I/I | PHAT | +5 |
| 4 | V/V | PHAT | +4 |
| 5 | gap/S | the gap lies between V (transmembrane) and G (outside), and S is outside: general open 3 | −3 |
| 6 | G/G | BLOSUM62 | +6 |

The higher penalties inside helices push gaps into the loops, where insertions and deletions
happen in real membrane proteins.

## Step 4: objective 2, lDDT

The second objective measures whether the alignment is **consistent with the structures**: if two
residues are aligned, their neighbourhoods should look the same in the two structures. It is a
reference-free, superposition-free version of the lDDT score (as in FoldMason):

1. For every ordered pair of proteins (query, target) and every column where both have a residue
   with structure, take the residue of the query as the *anchor*, and its neighbours closer than
   **R0 = 15 Å** in the query structure.
2. Every neighbour that the alignment maps to a residue of the target is scored by the difference
   between its distance to the anchor in the query and the distance between the corresponding
   residues in the target: 1 below 0.5 Å, 0.75 below 1 Å, 0.5 below 2 Å, 0.25 below 4 Å, and 0
   otherwise. Each neighbour is weighted by `1 − d / R0`, so close neighbours count more.
   Neighbours that the alignment does not map count in the weights only.
3. The anchor score is the weighted mean of its neighbours.
4. The anchor scores are combined into one value (`--lddt-aggregation`):
   - `HARMONIC_OVER_LENGTH` (default, the published definition): for every column, the harmonic
     mean of its mean anchor score and its occupancy (residues / N), summed over the columns and
     divided by the length of the alignment. It rewards compact alignments.
   - `PAIR_MEAN`: the mean of all anchor scores, independent of the length.

The value is between 0 and 1, and higher is better. Proteins without a usable structure (see the
mapping status in `metadata.json`) do not contribute.

In the example, moving the whole of `p2` one column to the right worsens both objectives:

```text
== 3. Objectives ==
seed alignment:
    PQTLLVAVLGNKE
    P-TLLIAVLG-KE
    PQSTLVAILGN-E
    SOPwTP 129.0   lDDT 0.9144   broken helices 0
p2 shifted one column:
    PQTLLVAVLGNKE
    P--TLLIAVLGKE
    PQSTLVAILGN-E
    SOPwTP 80.0   lDDT 0.8703   broken helices 0
```

The objectives often disagree on real data: the sequence score prefers some alignments and the
structures others. That is why the problem is solved as a multiobjective one, and the result is a
set of compromises instead of a single alignment.

> **Caution.** lDDT is almost blind to a helix aligned one or two positions out of register, and
> the default aggregation favours short alignments. Do not judge alignments by the second
> objective alone: compare them with a reference alignment as in
> [T4](T4-choosing-and-evaluating.md).

## Step 5: evaluation and gap-only columns

Before computing the objectives, every solution is **compacted**: the columns in which every row
has a gap are removed. Such columns appear after recombination and mutation, do not change the
alignment, and would only lengthen it.

```text
== 4. Gap-only columns ==
before evaluation:
    PQTLLVAVLG-NKE
    P-TLLIAVLG--KE
    PQSTLVAILG-N-E
after evaluation:
    PQTLLVAVLGNKE
    P-TLLIAVLG-KE
    PQSTLVAILGN-E
```

Both objectives are maximised. jMetal minimises, so the solutions store them negated; the output
files report them as computed. The evaluation of different solutions can run in parallel
(`--threads`) and does not change the results.

## Step 6: constraints

The problem has **no formal constraints** (in jMetal terms, `numberOfConstraints() = 0`): every
decoded alignment is valid. Not breaking transmembrane helices is a property of the operators
instead:

- a **broken helix** is a gap run between two consecutive transmembrane residues of a protein (the
  `broken helices` count of the example output);
- the crossovers never create one, and all mutations except the original `GAP_SHIFT` never create
  one either;
- the seed alignments can contain broken helices, and the search keeps them unless a mutation
  happens to move those gaps.

## Step 7: the initial population

The population starts with the **seed alignments** given with `--seeds`, in that order and without
duplicates. The rest is filled according to `--initialization`:

| Value | Members besides the seeds |
| --- | --- |
| `CROSSOVER` (default) | Children of the single-point crossover applied to random pairs of members |
| `CROSSOVER_DISTINCT` | The same, but a child equal to an existing member after compaction is discarded; the original fills 15–28 % of the population with such duplicates |
| `PERTURBED` | Crossover children followed by 10 mutations |

With MOEA/D the seeds are spread over equally spaced subproblems, so that they do not all start
at one end of the front.

## Step 8: crossover

The crossovers combine two parent alignments **without ever cutting a transmembrane helix**.

**Single point** (`--crossover SINGLE_POINT`, the default, class `TmSafeSinglePointCrossover`):

1. Choose a cut column `c` of the first parent at random among the **valid cuts**: columns that are
   not gap-only and where, for every protein, the residues on the left and the residue on the
   right are not two consecutive transmembrane residues.
2. The child takes columns `0..c` of the first parent.
3. For every protein, the child takes the rest of its residues **as aligned in the second parent**.
   The rest starts at different columns in the second parent for different proteins; proteins
   whose rest starts later are padded with gaps, so that the columns of the second parent stay
   aligned with each other.

Each crossover produces two children, (first, second) and (second, first). If no valid cut exists,
the child is a copy of its first parent.

```text
== 5. Crossover ==
first parent:
    PQTLLVAVLGNKE
    P-TLLIAVLG-KE
    PQSTLVAILGN-E
second parent:
    PQTLLVAVLGNKE--
    --PTLLIAVLGKE--
    PQSTLVAILGNE---
valid cuts of the first parent: [1, 2, 9, 10, 11]
single-point child, cut at column 1:
    PQTLLVAVLGNKE
    P--TLLIAVLGKE
    PQSTLVAILGNE-
    SOPwTP 72.0   lDDT 0.8717   broken helices 0
two-point child, cuts at columns 1 and 10:
    PQTLLVAVLGNKE
    P--TLLIAVLGKE
    PQSTLVAILGN-E
    SOPwTP 80.0   lDDT 0.8703   broken helices 0
```

Columns 3 to 8 are not valid cuts: they would separate residues of the helices. Column 12 is not
either, because the cut must leave at least one column on each side. In the single-point child,
`p2` keeps `P-` from the first parent and then needs its remaining residues `TLLIAVLGKE` from the
second parent, where they start at column 3 instead of column 2 as for the other proteins: one gap
of padding keeps them in register. The trailing columns that are gaps in every row disappear when
the child is evaluated.

**Two points** (`--crossover TWO_POINT`, class `TmSafeMultiPointCrossover`) repeats the operation
alternating the donor: the beginning comes from the first parent, the middle from the second, and
the end from the first again. Every cut is chosen among the valid cuts of the partial child after
the previous cut. With one point it gives exactly the same children as the single-point crossover.

## Step 9: mutation

A mutation changes one alignment. All mutations work on the gap layouts and keep every row a
valid encoding of its sequence. They are selected with `--mutation`; the probability of mutating
each child is `--mutation-probability`.

The output below applies three mutations to the same parent four times each (the random generator
is reset before every mutation):

```text
== 6. Mutations ==
parent:
    PQT-LLVAVLGNKE
    P-TLLIAVLG-KE-
    PQSTLVAILGN--E
    SOPwTP 36.0   lDDT 0.7941   broken helices 0
GAP_SHIFT (original), 4 children:
    PQT-LLVAVLGNKE
    P-TLLIAVLG-KE-
    PQSTLVAILGNE--
    SOPwTP 36.0   lDDT 0.7517   broken helices 0
    PQT-LLVAVLGNKE
    P-TLLIAVLG-KE-
    PQSTLVAILGNE--
    SOPwTP 36.0   lDDT 0.7517   broken helices 0
    PQT-LLVAVLGNKE
    P-TLLIAVLG-KE-
    P--QSTLVAILGNE
    SOPwTP -11.0   lDDT 0.7531   broken helices 0
    PQT-LLVAVLGNKE
    P-TLLIAVLG-KE-
    PQSTLVAILGNE--
    SOPwTP 36.0   lDDT 0.7517   broken helices 0
TM_SAFE_GAP_SHIFT, 4 children:
    PQT-LLVAVLGNKE
    PTLLIAVLG-KE--
    PQSTLVAILGN--E
    SOPwTP -6.0   lDDT 0.7461   broken helices 0
    PQT-LLVAVLGNKE
    P--TLLIAVLGKE-
    PQSTLVAILGN--E
    SOPwTP 52.0   lDDT 0.7719   broken helices 0
    PQT-LLVAVLGNKE
    P-TLLIAVLG-KE-
    PQSTLVAILGNE--
    SOPwTP 36.0   lDDT 0.7517   broken helices 0
    PQT-LLVAVLGNKE
    P-TLLIAVLG--KE
    PQSTLVAILGN--E
    SOPwTP 56.0   lDDT 0.7453   broken helices 0
MERGE_GAP_RUNS, 4 children:
    PQT-LLVAVLGNKE
    P--TLLIAVLGKE-
    PQSTLVAILGN--E
    SOPwTP 52.0   lDDT 0.7719   broken helices 0
    PQT-LLVAVLGNKE
    P-TLLIAVLG--KE
    PQSTLVAILGN--E
    SOPwTP 56.0   lDDT 0.7453   broken helices 0
    PQT-LLVAVLGNKE
    P-TLLIAVLGKE--
    PQSTLVAILGN--E
    SOPwTP 44.0   lDDT 0.7335   broken helices 0
    PQT-LLVAVLGNKE
    P-TLLIAVLG--KE
    PQSTLVAILGN--E
    SOPwTP 56.0   lDDT 0.7453   broken helices 0
```

**`GAP_SHIFT`** (default, class `TmAwareGapShiftMutation`) is the mutation of the published
experiments:

1. Choose at random a protein that has a run of at least two gaps, and one of its runs.
2. Remove the run and propose a new position at random. Accept it with a probability that depends
   on what occupies that position in the row: 0.60 a gap, 0.95 a residue without structure, 0.05 a
   transmembrane residue, 0.90 any other residue. After 60 rejections, a position chosen uniformly
   at random is used.
3. Insert the run at the accepted position; the length of the alignment does not change.

Only `p3` has a run of two gaps in the parent, so it is the only protein that this mutation
changes. Helices are only discouraged, not protected: a run can still land inside a helix, or a gap
run already inside a helix can grow.

**`TM_SAFE_GAP_SHIFT`** (class `TmSafeGapShiftMutation`) moves a gap run in the same way, with
three differences: the new position is only chosen among the positions that **do not separate two
transmembrane residues**, any protein with gaps can be chosen (`p2` and `p3` change in the
example), and the boundary is judged by the residues on both sides, as in the crossover. The
acceptance probabilities among valid positions are kept.

**`MERGE_GAP_RUNS`** (class `MergeAdjacentGapRunsMutation`) joins two consecutive gap runs of one
protein: one run moves next to the other, shifting the residues between them. It creates no new
gap boundary, so it never breaks a helix; when one of the two runs breaks a helix, that is the one
moved, which repairs the break. The length of the alignment does not change. In the first child,
the run of `p2` at column 10 moves to column 2.

**`RANDOM_CHOICE`** (class `RandomChoiceMutation`) applies `TM_SAFE_GAP_SHIFT` or
`MERGE_GAP_RUNS`, chosen uniformly at random each time the child is mutated.

**How much a mutation changes an alignment** matters as much as what it does. Applied once to each
seed alignment of `classF`:

| Mutation | Children dominated by their parent | Mean change of SOPwTP | Mean change of lDDT |
| --- | --- | --- | --- |
| `TM_SAFE_GAP_SHIFT` | 82 % | −8 378 | −0.0047 |
| `MERGE_GAP_RUNS` | 50 % | −687 | −0.0003 |
| Insert a gap in every row (removed) | 100 % | −44 185 | −0.0331 |

The last row is a mutation of TM-MSAligner that inserted one gap in every row, each at an
independent position. Between the gaps of two rows, the residues of one row are shifted with
respect to the other: on average a third of the alignment for every pair of proteins. Its children
never survived the selection, so it did not take part in the search, and it was removed together
with its variant that opened the gap between two residues.

## Step 10: reproducibility

Every random decision (initial population, choice of parents, cuts, mutation positions) is drawn
from one generator seeded with `--seed`, through the interface `operator.RandomSource`. The same
version of the tool with the same options and seed produces the same result, whatever the number
of threads. The tests replace the generator with fixed sequences to check every operator step by
step.

## What next

- [T5](T5-comparing-algorithms-and-operators.md) compares algorithms and operators on real data.
- [T7](T7-adding-an-operator.md) shows how to add a new operator.
- The [command-line reference](../reference/cli.md) lists every option.
