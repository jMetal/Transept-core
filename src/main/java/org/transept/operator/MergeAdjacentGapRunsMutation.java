package org.transept.operator;

import java.util.ArrayList;
import java.util.List;
import org.transept.alignment.GapIntervals;
import org.transept.alignment.MsaSolution;
import org.transept.model.ProteinSequence;
import org.uma.jmetal.operator.mutation.MutationOperator;

/**
 * Merges two consecutive gap runs of one row into a single run.
 *
 * <p>TM-safe version of {@code MergeAdjunctedGapsGroupsMSAMutation} (MAGG) of M2Align and
 * TM-MSAligner. With the mutation probability, a row with at least two gap runs is chosen at
 * random, then a pair of consecutive runs. One run of the pair is moved next to the other, so
 * the residues between them shift and the merged run sits at the boundary of the run that stays.
 * No new boundary is created, so the operator never breaks a helix; if exactly one of the two
 * boundaries lies inside a helix, the run there is the one moved, which repairs the break. When
 * both or neither boundary is valid, the run that stays is chosen at random. The alignment keeps
 * its length.
 *
 * <p>Unlike the original, which looped forever when no row had two runs, the solution is left
 * unchanged in that case.
 */
public final class MergeAdjacentGapRunsMutation implements MutationOperator<MsaSolution> {

  private final double probability;
  private final RandomSource random;

  /**
   * Creates the operator using the global jMetal generator.
   *
   * @param probability probability of mutating a solution
   */
  public MergeAdjacentGapRunsMutation(double probability) {
    this(probability, RandomSource.jmetal());
  }

  /**
   * Creates the operator.
   *
   * @param probability probability of mutating a solution
   * @param random source of random numbers
   */
  public MergeAdjacentGapRunsMutation(double probability, RandomSource random) {
    if (!(probability >= 0.0 && probability <= 1.0)) {
      throw new IllegalArgumentException("Mutation probability must be in [0, 1]: " + probability);
    }
    this.probability = probability;
    this.random = random;
  }

  @Override
  public MsaSolution execute(MsaSolution solution) {
    if (solution == null) {
      throw new IllegalArgumentException("Solution cannot be null");
    }
    if (random.nextDouble() < probability) {
      List<Integer> eligible = new ArrayList<>();
      for (int sequence = 0; sequence < solution.variables().size(); sequence++) {
        if (solution.variables().get(sequence).runCount() > 1) {
          eligible.add(sequence);
        }
      }
      if (!eligible.isEmpty()) {
        int sequence = eligible.get(random.nextInt(0, eligible.size() - 1));
        GapIntervals row = solution.variables().get(sequence);
        int left = random.nextInt(0, row.runCount() - 2);
        solution
            .variables()
            .set(sequence, merge(solution.proteins().get(sequence), row, left));
      }
    }
    return solution;
  }

  @Override
  public double mutationProbability() {
    return probability;
  }

  private GapIntervals merge(ProteinSequence protein, GapIntervals row, int left) {
    boolean[] gaps = row.toMask();
    int[] residuesBefore = TmBoundaries.residuesBefore(gaps);
    boolean leftValid = !TmBoundaries.splitsHelix(protein, residuesBefore[row.start(left)]);
    boolean rightValid = !TmBoundaries.splitsHelix(protein, residuesBefore[row.start(left + 1)]);
    boolean keepLeft = leftValid == rightValid ? random.nextDouble() < 0.5 : leftValid;
    boolean[] merged;
    if (keepLeft) {
      int length = row.runLength(left + 1);
      merged = TmBoundaries.removeColumns(gaps, row.start(left + 1), length);
      merged = TmBoundaries.insertGaps(merged, row.end(left) + 1, length);
    } else {
      int length = row.runLength(left);
      merged = TmBoundaries.removeColumns(gaps, row.start(left), length);
      merged = TmBoundaries.insertGaps(merged, row.start(left + 1) - length, length);
    }
    return GapIntervals.fromMask(merged);
  }
}
