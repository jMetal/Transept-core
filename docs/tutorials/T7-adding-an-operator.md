# T7. Adding a Variation Operator

| | |
| --- | --- |
| **Level** | Advanced (developers) |
| **Version** | 1.0 (2026-11-04) |
| **Time** | about 30 minutes to read; implementing and evaluating an operator takes longer |
| **Prerequisites** | [T2. Problem formulation](T2-problem-formulation.md), Java 21, the coding guidelines of the project (`JAVA_CODING_GUIDELINES.md`) |

This tutorial explains how to add a mutation or a crossover to TM-M2StructAlign: the interfaces
it implements, the rules every operator must follow, the tests it needs, how to make it available
in the command line, and how to decide whether it is worth keeping. The worked example is
`MergeAdjacentGapRunsMutation`, a short operator of the project.

## Step 1: the interfaces

Operators are jMetal operators on `MsaSolution`:

```java
public interface MutationOperator<S> {
  S execute(S solution);
  double mutationProbability();
}

public interface CrossoverOperator<S> {
  List<S> execute(List<S> parents);
  double crossoverProbability();
  int numberOfRequiredParents();
  int numberOfGeneratedChildren();
}
```

An `MsaSolution` holds one `GapIntervals` per protein (`solution.variables()`), the proteins with
their topologies (`solution.proteins()`), and the objective values. A mutation changes the solution
it receives and returns it; a crossover returns new solutions and leaves its parents unchanged.

## Step 2: the rules

Every operator of the project follows five rules. The tests of Step 4 check them.

1. **Valid encodings only.** `GapIntervals` is immutable: build a new layout and store it with
   `solution.variables().set(sequence, layout)`. Every row must keep its residues and all rows
   the same length. Build layouts with `GapIntervals.fromMask(boolean[])` or
   `GapIntervals.of(length, start, end, ...)`, which validate them.
2. **Never break a helix.** A gap run must not end up between two consecutive transmembrane
   residues. The package-private class `operator.TmBoundaries` has the helpers:
   `residuesBefore(mask)` counts the residues before every insertion position, and
   `splitsHelix(protein, residues)` tells whether gaps placed there would break a helix. Gaps inside
   helices inherited from the seeds may stay; the operator must not create new ones.
3. **All randomness from `RandomSource`.** Take a `RandomSource` in the constructor and draw every
   random number from it; offer a convenience constructor with `RandomSource.jmetal()`. Never use
   `Math.random()` or `new Random()`: runs would not be reproducible from `--seed`, and the tests
   could not control the operator.
4. **Do nothing gracefully.** If the operator cannot act (no gap run to move, no valid position),
   return the solution unchanged. Never loop until a condition holds: two operators of the original
   code looped forever on alignments where the condition could not hold.
5. **No shared mutable state.** The evaluation runs in parallel; operators run in one thread, but
   should not keep state between calls besides the random source.

## Step 3: the example

`MergeAdjacentGapRunsMutation` joins two consecutive gap runs of one row. Its core, simplified:

```java
@Override
public MsaSolution execute(MsaSolution solution) {
  if (random.nextDouble() < probability) {
    List<Integer> eligible = rowsWithAtLeastTwoRuns(solution);       // rule 4: may be empty
    if (!eligible.isEmpty()) {
      int sequence = eligible.get(random.nextInt(0, eligible.size() - 1));
      GapIntervals row = solution.variables().get(sequence);
      int left = random.nextInt(0, row.runCount() - 2);
      solution.variables().set(sequence,                               // rule 1: new layout
          merge(solution.proteins().get(sequence), row, left));
    }
  }
  return solution;
}

private GapIntervals merge(ProteinSequence protein, GapIntervals row, int left) {
  boolean[] gaps = row.toMask();
  int[] residuesBefore = TmBoundaries.residuesBefore(gaps);
  boolean leftValid = !TmBoundaries.splitsHelix(protein, residuesBefore[row.start(left)]);
  boolean rightValid = !TmBoundaries.splitsHelix(protein, residuesBefore[row.start(left + 1)]);
  boolean keepLeft = leftValid == rightValid ? random.nextDouble() < 0.5 : leftValid;  // rule 2
  ...  // move one run next to the other with TmBoundaries.removeColumns and insertGaps
  return GapIntervals.fromMask(merged);
}
```

The operator creates no new gap boundary, so it cannot break a helix; when one of the two runs
breaks one, that is the run it moves.

## Step 4: the tests

Each operator has a unit test class in `src/test/java/org/tmm2structalign/operator/`, following
the conventions of the project (JUnit Jupiter, `given_when_then` names, `@DisplayName`,
Arrange/Act/Assert). The utilities in `org.tmm2structalign.testing` help:

| Utility | Use |
| --- | --- |
| `Alignments.proteins("ACDEF:OMMMO", ...)` | Proteins from residues and topology |
| `Alignments.solution(proteins, "AC-DEF", ...)` | A solution from aligned rows |
| `Alignments.rows(solution)` | The aligned rows of a solution, to compare with strings |
| `Alignments.randomSolution(proteins, maxGaps, random)` | Random valid alignments, often with broken helices |
| `Alignments.helixBreakingRuns(protein, row)` | Number of gap runs that break a helix |
| `SeededRandom(seed)` | A reproducible `RandomSource` |

Write at least these tests:

1. **Exact behaviour with a mocked random source.** Mock `RandomSource` with Mockito, fix every draw,
   and compare the result with the expected rows. One test per case of the operator (each
   direction, each fallback). This is the only way to check the operator step by step.
2. **Property test on random alignments.** Apply the operator hundreds of times to random
   alignments of proteins with helices, and check that the number of helix-breaking runs never
   grows, the rows keep their residues (`solution.toMsa()` validates them), and the length behaves
   as documented.
3. **Edge cases.** The operator cannot act (rule 4), probability 0, invalid arguments.

The integration test `TmSafetyIT` runs NSGA-II with every TM-safe mutation and crossover of the
command line on a synthetic dataset and checks that no final alignment breaks a helix. It covers
every value of `MutationType` except `GAP_SHIFT`, so a new mutation is tested automatically; keep
`MutationType.isTmSafe()` consistent with it.

## Step 5: make it available

1. Add a value to `algorithm.MutationType` (or `CrossoverType`), with a Javadoc comment that says
   what the operator does.
2. Build the operator in `AlgorithmFactory.mutation(...)` (or `crossover()`). The `switch` is
   exhaustive, so the compiler points to every place to update.
3. Update the description of `--mutation` in `cli.RunCommand`.
4. Document it: the table of `docs/reference/cli.md`, the operator section of
   [T2](T2-problem-formulation.md), and, if it helps to show it, an entry in
   `example.FormulationTour` (then paste the new output into T2: `FormulationTourTest` fails
   until you do).
5. Keep the defaults unchanged: they determine published results, and changing them is a decision
   of the authors.

## Step 6: decide whether it is worth keeping

An operator that is correct can still be useless. Evaluate it in two stages before proposing it.

**Stage 1: how much it changes an alignment.** Apply it once to each seed alignment of a dataset,
many times, evaluate the children, and count how many are dominated by their parent and how much
each objective changes on average. A mutation whose children are almost always much worse never
survives the selection, and does not take part in the search. This is what happened to a mutation
of the original code that inserted a gap in every row at independent positions: on `classF`, all
its children were dominated, with a mean loss of about 44 000 of SOPwTP, against about 8 000 for
`TM_SAFE_GAP_SHIFT`.

**Stage 2: whether it helps the search.** Run NSGA-II with and without the operator, with several
random seeds (at least 6, better 10 or more) on two or more datasets, and compare the
hypervolume of the results, normalised over all the runs of a dataset, with a paired test such as
Wilcoxon's signed-rank test. A useful way to add a mutation without removing the existing one is
`RandomChoiceMutation`, which applies one of several mutations chosen at random. `scripts/pareto_analysis.py`
computes the hypervolume of every run ([T5](T5-comparing-algorithms-and-operators.md)).

A low impact in stage 1 is necessary but not sufficient. A mutation that inserted gaps in two groups
of rows a few columns apart reduced the mean loss of SOPwTP to about 250, but added to
`TM_SAFE_GAP_SHIFT` it lowered the hypervolume on `classF` in all six seeds of the comparison, and it
was removed. The same happened to a mutation that moved a whole helix a few columns.

## Step 7: a crossover

A crossover follows the same rules. Reuse `TmSafeSinglePointCrossover`: its public methods
`compatibleCuts(first)` (the cuts that do not split a helix) and `recombine(first, second, cut)`
(the child of one cut) are the building blocks of `TmSafeMultiPointCrossover`, which applies
several cuts in sequence. A crossover must return copies when it does not recombine, and must check
that both parents align the same protein set.
