# T2. Problem Formulation, Encoding and Variation Operators

| | |
| --- | --- |
| **Level** | Introductory |
| **Version** | 1.0 (2026-11-04) |
| **Time** | about 45 minutes |
| **Prerequisites** | [T1. First run in five minutes](T1-first-run.md) |
| **Companion** | `org.transept.example.FormulationTour` |

This tutorial describes what Transept optimises and how it does it: the inputs, the two
objectives, the constraints, how an alignment is encoded, how the initial population is built,
and what every crossover and mutation does. A small example with three proteins runs through it.
Its output is printed by a class of the project, so you can run it and change it:

```bash
java -cp target/transept-core-0.1.0-jar-with-dependencies.jar \
  org.transept.example.FormulationTour
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
4. The objective is the mean of all anchor scores. It does not depend on the length of the
   alignment.

The original implementation combined the anchor scores of every column with the column occupancy
(their harmonic mean), summed over the columns and divided by the alignment length. That rewarded
compact alignments: on the precomputed alignments of the benchmark it was *negatively* correlated
with the agreement with the GPCRdb references in 4 of 5 datasets, and it scored the misaligned
outputs of Caretta among the best. The mean of the anchor scores is positively correlated with
the references in all of them, and gives better alignments when it guides the search.

The value is between 0 and 1, and higher is better. Proteins without a usable structure (see the
mapping status in `metadata.json`) do not contribute.

In the example, moving the whole of `p2` one column to the right worsens both objectives:

```text
== 3. Objectives ==
seed alignment:
    PQTLLVAVLGNKE
    P-TLLIAVLG-KE
    PQSTLVAILGN-E
    SOPwTP 129.0   lDDT 0.9167   broken helices 0
p2 shifted one column:
    PQTLLVAVLGNKE
    P--TLLIAVLGKE
    PQSTLVAILGN-E
    SOPwTP 80.0   lDDT 0.8393   broken helices 0
```

The objectives often disagree on real data: the sequence score prefers some alignments and the
structures others. That is why the problem is solved as a multiobjective one, and the result is a
set of compromises instead of a single alignment.

> **Caution.** A superposition-free lDDT is almost blind to a helix aligned one or two positions
> out of register: the distances along a helix repeat every turn. Do not judge alignments by the
> second objective alone: compare them with a reference alignment as in
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
- neither the crossover nor the mutations ever create one (the mutation of the original
  implementation only discouraged them);
- the seed alignments can contain broken helices, and the search keeps them unless a mutation
  happens to move those gaps.

## Step 7: the initial population

The population starts with the **seed alignments** given with `--seeds`, in that order and without
duplicates. The rest is filled with children of the single-point crossover (Step 8) applied to
random pairs of members. A child equal to an existing member once its gap-only columns are removed
is discarded: without this check, 15–28 % of the population were such duplicates. If the
crossover stops producing new alignments (few, very similar seeds), the remaining members are
copies of random members changed by 10 mutations.

With MOEA/D the seeds are spread over equally spaced subproblems, so that they do not all start
at one end of the front.

## Step 8: crossover

The crossover (class `TmSafeSinglePointCrossover`) combines two parent alignments **without ever
cutting a transmembrane helix**:

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
    SOPwTP 72.0   lDDT 0.8410   broken helices 0
```

Columns 3 to 8 are not valid cuts: they would separate residues of the helices. Column 12 is not
either, because the cut must leave at least one column on each side. In the child,
`p2` keeps `P-` from the first parent and then needs its remaining residues `TLLIAVLGKE` from the
second parent, where they start at column 3 instead of column 2 as for the other proteins: one gap
of padding keeps them in register. The trailing columns that are gaps in every row disappear when
the child is evaluated.

**How much padding a cut needs varies a lot.** With few gaps it can be zero, but in alignments with
long variable regions a random cut can displace a block of the second parent by hundreds of columns,
which also makes alignments grow. With `--cut-candidates k`, the crossover draws `k` valid cuts and
uses the one that adds the fewest gaps. The default, 1, is a cut at random. Preferring those cuts
does not always help, because it also explores less: compare both on your data.

A two-point variant was also evaluated; it did not improve the results and was removed.

**`TM_SAFE_UNIFORM`** (class `TmSafeUniformCrossover`, `--crossover TM_SAFE_UNIFORM`) is the uniform
crossover of SAGA (Notredame and Higgins, 1996). It looks for the **consistent positions** of the two
parents: places where, in every protein, both parents have the same residues on the left. Between
two consecutive consistent positions, each parent holds the same residues of every protein, aligned
in its own way, so a child can take that block from either parent:

```text
first parent      second parent     first child       second child
A-CD-E            AC-DE-            A-CDE-            AC-D-E
AC-DE-            A-CD-E            AC-D-E            A-CDE-
```

Here the consistent positions are after columns 1, 3 and 4 and at both ends, which makes four
blocks: `A`, the next two columns, `D` and the last two columns; the first and third are the same in
both parents. The first child takes the second block from the first parent and the rest from the
second; the second child does the opposite. No gap is added: every block keeps the gaps it had in its
parent, although a child can be longer than both parents if it takes the longer version of several
blocks. A child separates two transmembrane residues only where one of its parents did. When the parents agree only at the
start and the end, the children are copies of them. One probability draw decides whether the two
parents are crossed. The initial population is always filled with the single-point crossover.

## Step 9: mutation

A mutation changes one alignment. All mutations work on the gap layouts and keep every row a
valid encoding of its sequence. They are selected with `--mutation`; the probability of mutating
each child is `--mutation-probability`.

The output below applies the two mutations to the same parent four times each (the random generator
is reset before every mutation):

```text
== 6. Mutations ==
parent:
    PQT-LLVAVLGNKE
    P-TLLIAVLG-KE-
    PQSTLVAILGN--E
    SOPwTP 36.0   lDDT 0.7543   broken helices 0
TM_SAFE_GAP_SHIFT, 4 children:
    PQT-LLVAVLGNKE
    PTLLIAVLG-KE--
    PQSTLVAILGN--E
    SOPwTP -6.0   lDDT 0.7410   broken helices 0
    PQT-LLVAVLGNKE
    P--TLLIAVLGKE-
    PQSTLVAILGN--E
    SOPwTP 52.0   lDDT 0.7357   broken helices 0
    PQT-LLVAVLGNKE
    P-TLLIAVLG-KE-
    PQSTLVAILGNE--
    SOPwTP 36.0   lDDT 0.7482   broken helices 0
    PQT-LLVAVLGNKE
    P-TLLIAVLG--KE
    PQSTLVAILGN--E
    SOPwTP 56.0   lDDT 0.7399   broken helices 0
MERGE_GAP_RUNS, 4 children:
    PQT-LLVAVLGNKE
    P--TLLIAVLGKE-
    PQSTLVAILGN--E
    SOPwTP 52.0   lDDT 0.7357   broken helices 0
    PQT-LLVAVLGNKE
    P-TLLIAVLG--KE
    PQSTLVAILGN--E
    SOPwTP 56.0   lDDT 0.7399   broken helices 0
    PQT-LLVAVLGNKE
    P-TLLIAVLGKE--
    PQSTLVAILGN--E
    SOPwTP 44.0   lDDT 0.7268   broken helices 0
    PQT-LLVAVLGNKE
    P-TLLIAVLG--KE
    PQSTLVAILGN--E
    SOPwTP 56.0   lDDT 0.7399   broken helices 0
```

**`TM_SAFE_GAP_SHIFT`** (default, class `TmSafeGapShiftMutation`) moves a run of gaps of one
protein:

1. Choose at random a protein that has gaps, and one of its gap runs.
2. Remove the run and propose a new position at random among the positions that **do not
   separate two transmembrane residues** (judged by the residues on both sides, as in the
   crossover). Accept it with a probability that depends on what occupies that position in the
   row: 0.60 a gap, 0.95 a residue without structure, 0.90 any other residue. After 60
   rejections, a valid position chosen uniformly at random is used.
3. Insert the run at the accepted position; the length of the alignment does not change.

By default the new position can be anywhere in the row, so a run can jump over many residues, even a
whole helix. `--max-shift n` limits a move to `n` columns.

In the example, `p2` and `p3` change. The original implementation used a variant that only
*discouraged* positions inside helices (acceptance 0.05), so a run could still break a helix; it
was replaced by this one, which gives alignments of the same quality without ever breaking one.

**`MERGE_GAP_RUNS`** (class `MergeAdjacentGapRunsMutation`) joins two consecutive gap runs of one
protein: one run moves next to the other, shifting the residues between them. It creates no new
gap boundary, so it never breaks a helix; when one of the two runs breaks a helix, that is the one
moved, which repairs the break. The length of the alignment does not change. In the first child,
the run of `p2` at column 10 moves to column 2.

**`RANDOM_CHOICE`** (class `RandomChoiceMutation`) applies `TM_SAFE_GAP_SHIFT` or
`MERGE_GAP_RUNS`, chosen uniformly at random each time the child is mutated.

**`TM_SAFE_BLOCK_SHIFT`** (class `TmSafeBlockShiftMutation`) moves a vertical block of gaps in
several proteins at once, based on the block shuffling of SAGA (Notredame and Higgins, 1996):

1. Choose at random a protein that has gaps, and one of its gap runs; its columns define the block.
2. The other proteins with gaps in all those columns move with it: all of them half of the time,
   otherwise a random subset, so that a block can also be split.
3. Move the block left or right by up to `blockShift` columns (1 by default). A protein where the
   new position would separate two transmembrane residues stays where it is. The length of the
   alignment does not change.

Unlike `TM_SAFE_GAP_SHIFT`, the proteins that move keep their alignment with each other, so the
change is small and local.

**`TM_SAFE_WINDOW_REALIGN`** (class `TmSafeWindowRealignMutation`) realigns one protein inside a
window of columns, in the spirit of the local realignment of iterative aligners:

1. Choose at random a protein that has gaps, one of its gap runs, and a window of `realignWindow`
   columns (50 by default) that contains the run.
2. Place the residues of that protein inside the window again, by dynamic programming, against the
   other proteins, keeping their order: the score is the substitution score (PHAT inside the membrane, BLOSUM62 outside)
   and the gap penalties of SOPwTP. A placement that would separate two transmembrane residues is
   not allowed. The structures are not used.
3. The columns outside the window, the other proteins and the length of the alignment do not change.

`blockShift` and `realignWindow` are set in a configuration ([configuration](../reference/configuration.md));
with `--mutation`, the defaults are used. `RANDOM_CHOICE` does not include these two mutations.

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
