package org.tmm2structalign.operator;

import java.util.ArrayList;
import java.util.List;
import org.tmm2structalign.alignment.GapIntervals;
import org.tmm2structalign.alignment.MsaSolution;
import org.tmm2structalign.model.ProteinSequence;
import org.uma.jmetal.operator.mutation.MutationOperator;

/**
 * Inserts gaps into every row at nearby positions, so that only a short stretch of the
 * alignment is shifted, without breaking a transmembrane helix.
 *
 * <p>With the mutation probability, an insertion position {@code c} of the alignment and an
 * offset {@code d}, with {@code 1 <= |d| <= window}, are chosen at random, and the rows are
 * split at random into two non-empty groups. Every row of the first group receives
 * {@code gapLength} gaps at {@code c}; every row of the second group, at {@code c + d}. Between
 * the two positions the residues of one group are shifted with respect to the other, by at most
 * {@code window} columns, whereas {@link TmSafeInsertGapMutation}, which places the gaps of
 * every row independently, shifts about a third of the alignment for each pair of rows. If a
 * target position would break a helix, the nearest valid position is used (the end of a row is
 * always valid). It is the gap insertion of SAGA (Notredame and Higgins, 1996) with random
 * groups instead of a guide tree. The alignment grows by {@code gapLength} columns.
 */
public final class TmSafeLocalGapInsertionMutation implements MutationOperator<MsaSolution> {

  /** Default maximum distance between the insertion positions of the two groups. */
  public static final int DEFAULT_WINDOW = 3;

  private final double probability;
  private final int gapLength;
  private final int window;
  private final RandomSource random;

  /**
   * Creates the operator with one gap per row, the default window and the global jMetal
   * generator.
   *
   * @param probability probability of mutating a solution
   */
  public TmSafeLocalGapInsertionMutation(double probability) {
    this(probability, 1, DEFAULT_WINDOW, RandomSource.jmetal());
  }

  /**
   * Creates the operator.
   *
   * @param probability probability of mutating a solution
   * @param gapLength gaps inserted into every row
   * @param window maximum distance between the insertion positions of the two groups
   * @param random source of random numbers
   */
  public TmSafeLocalGapInsertionMutation(
      double probability, int gapLength, int window, RandomSource random) {
    if (!(probability >= 0.0 && probability <= 1.0)) {
      throw new IllegalArgumentException("Mutation probability must be in [0, 1]: " + probability);
    }
    if (gapLength < 1) {
      throw new IllegalArgumentException("gapLength must be positive: " + gapLength);
    }
    if (window < 1) {
      throw new IllegalArgumentException("window must be positive: " + window);
    }
    this.probability = probability;
    this.gapLength = gapLength;
    this.window = window;
    this.random = random;
  }

  @Override
  public MsaSolution execute(MsaSolution solution) {
    if (solution == null) {
      throw new IllegalArgumentException("Solution cannot be null");
    }
    int rows = solution.variables().size();
    if (rows > 1 && random.nextDouble() < probability) {
      int length = solution.alignmentLength();
      int first = random.nextInt(0, length);
      int offset = random.nextInt(1, window);
      int second = random.nextDouble() < 0.5 ? first - offset : first + offset;
      if (second < 0 || second > length) {
        second = Math.max(0, Math.min(length, 2 * first - second));
      }
      List<Boolean> inFirstGroup = randomGroups(rows);
      for (int sequence = 0; sequence < rows; sequence++) {
        int target = inFirstGroup.get(sequence) ? first : second;
        boolean[] gaps = solution.variables().get(sequence).toMask();
        int position = nearestValid(solution.proteins().get(sequence), gaps, target);
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

  /** Returns the maximum distance between the insertion positions of the two groups. */
  public int window() {
    return window;
  }

  /**
   * Returns the valid insertion position closest to the target; on a tie, the one on the left.
   *
   * @param protein protein of the row
   * @param gaps gap mask of the row
   * @param target preferred insertion position, in {@code [0, gaps.length]}
   * @return a position whose boundary does not split a helix
   */
  static int nearestValid(ProteinSequence protein, boolean[] gaps, int target) {
    int[] residuesBefore = TmBoundaries.residuesBefore(gaps);
    int position = gaps.length;
    boolean found = false;
    for (int distance = 0; !found && distance <= gaps.length; distance++) {
      for (int candidate : new int[] {target - distance, target + distance}) {
        if (!found
            && candidate >= 0
            && candidate <= gaps.length
            && !TmBoundaries.splitsHelix(protein, residuesBefore[candidate])) {
          position = candidate;
          found = true;
        }
      }
    }
    return position;
  }

  /** Splits the rows into two non-empty groups, each row joining the first with probability 1/2. */
  private List<Boolean> randomGroups(int rows) {
    List<Boolean> inFirstGroup = new ArrayList<>(rows);
    int members = 0;
    for (int sequence = 0; sequence < rows; sequence++) {
      boolean member = random.nextDouble() < 0.5;
      inFirstGroup.add(member);
      members += member ? 1 : 0;
    }
    if (members == 0 || members == rows) {
      int flipped = random.nextInt(0, rows - 1);
      inFirstGroup.set(flipped, !inFirstGroup.get(flipped));
    }
    return inFirstGroup;
  }
}
