package org.tmm2structalign.operator;

import java.util.ArrayList;
import java.util.List;
import org.tmm2structalign.alignment.GapIntervals;
import org.tmm2structalign.alignment.MsaSolution;
import org.tmm2structalign.model.ProteinSequence;
import org.uma.jmetal.operator.mutation.MutationOperator;

/**
 * Base of the mutations that insert the same number of gaps into every row, so that the rows
 * keep a common length, at positions that never break a transmembrane helix.
 *
 * <p>With the mutation probability, every row receives {@code gapLength} gaps at an insertion
 * position chosen uniformly among the <i>preferred</i> positions defined by the subclass that are
 * also valid ({@link TmBoundaries}). If a row has no such position, any valid position is used;
 * the end of the row is always valid. The alignment grows by {@code gapLength} columns.
 */
abstract class TmSafeGapInsertion implements MutationOperator<MsaSolution> {

  private final double probability;
  private final int gapLength;
  private final RandomSource random;

  TmSafeGapInsertion(double probability, int gapLength, RandomSource random) {
    if (!(probability >= 0.0 && probability <= 1.0)) {
      throw new IllegalArgumentException("Mutation probability must be in [0, 1]: " + probability);
    }
    if (gapLength < 1) {
      throw new IllegalArgumentException("gapLength must be positive: " + gapLength);
    }
    this.probability = probability;
    this.gapLength = gapLength;
    this.random = random;
  }

  /**
   * Tells whether an insertion position is preferred.
   *
   * @param gaps gap mask of the row
   * @param position insertion position, in {@code [0, gaps.length]}
   * @return true if the subclass prefers to insert there
   */
  abstract boolean isPreferred(boolean[] gaps, int position);

  @Override
  public MsaSolution execute(MsaSolution solution) {
    if (solution == null) {
      throw new IllegalArgumentException("Solution cannot be null");
    }
    if (random.nextDouble() < probability) {
      for (int sequence = 0; sequence < solution.variables().size(); sequence++) {
        boolean[] gaps = solution.variables().get(sequence).toMask();
        int position = choosePosition(solution.proteins().get(sequence), gaps);
        solution
            .variables()
            .set(
                sequence,
                GapIntervals.fromMask(TmBoundaries.insertGaps(gaps, position, gapLength)));
      }
    }
    return solution;
  }

  @Override
  public double mutationProbability() {
    return probability;
  }

  /** Returns the number of gaps inserted into every row. */
  public int gapLength() {
    return gapLength;
  }

  private int choosePosition(ProteinSequence protein, boolean[] gaps) {
    int[] residuesBefore = TmBoundaries.residuesBefore(gaps);
    List<Integer> valid = new ArrayList<>();
    List<Integer> preferred = new ArrayList<>();
    for (int position = 0; position <= gaps.length; position++) {
      if (!TmBoundaries.splitsHelix(protein, residuesBefore[position])) {
        valid.add(position);
        if (isPreferred(gaps, position)) {
          preferred.add(position);
        }
      }
    }
    List<Integer> candidates = preferred.isEmpty() ? valid : preferred;
    return candidates.get(random.nextInt(0, candidates.size() - 1));
  }
}
