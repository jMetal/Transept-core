package org.transept.operator;

import java.util.ArrayList;
import java.util.List;
import org.transept.alignment.GapIntervals;
import org.transept.alignment.MsaSolution;
import org.transept.model.ProteinSequence;
import org.uma.jmetal.operator.mutation.MutationOperator;

/**
 * Moves a vertical block of gaps a few columns, in several rows at once, without ever breaking a
 * transmembrane helix.
 *
 * <p>With the mutation probability, a gap run of a row is chosen at random. Its columns define the
 * block, and the rows with gaps in all those columns can move with it: all of them half of the
 * time, otherwise a random subset, so that a block can also be split. The block moves left or right
 * by a random number of columns up to {@code maxShift}. A row where the new position would separate
 * two transmembrane residues, or would leave the row, stays where it is. If no row can move in that
 * direction the other one is tried. The alignment keeps its length.
 *
 * <p>Moving the rows together keeps them aligned with each other, unlike moving a run of one row.
 * It is based on the block shuffling of SAGA (Notredame and Higgins, 1996), without its guide tree.
 */
public final class TmSafeBlockShiftMutation implements MutationOperator<MsaSolution> {

  /** Default largest number of columns a block can move, as in SAGA. */
  public static final int DEFAULT_MAX_SHIFT = 1;

  private final double probability;
  private final int maxShift;
  private final RandomSource random;

  /**
   * Creates the operator.
   *
   * @param probability probability of mutating a solution
   * @param maxShift largest number of columns a block can move
   * @param random source of random numbers
   */
  public TmSafeBlockShiftMutation(double probability, int maxShift, RandomSource random) {
    if (!(probability >= 0.0 && probability <= 1.0)) {
      throw new IllegalArgumentException("Mutation probability must be in [0, 1]: " + probability);
    }
    if (maxShift < 1) {
      throw new IllegalArgumentException("maxShift must be positive: " + maxShift);
    }
    this.probability = probability;
    this.maxShift = maxShift;
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
        if (solution.variables().get(sequence).runCount() > 0) {
          eligible.add(sequence);
        }
      }
      if (!eligible.isEmpty()) {
        int sequence = eligible.get(random.nextInt(0, eligible.size() - 1));
        GapIntervals row = solution.variables().get(sequence);
        int run = random.nextInt(0, row.runCount() - 1);
        moveBlock(solution, sequence, row.start(run), row.runLength(run));
      }
    }
    return solution;
  }

  @Override
  public double mutationProbability() {
    return probability;
  }

  private void moveBlock(MsaSolution solution, int chosen, int start, int size) {
    List<Integer> rows = rowsOfBlock(solution, chosen, start, size);
    int distance = random.nextInt(1, maxShift);
    int direction = random.nextDouble() < 0.5 ? -1 : 1;
    if (!move(solution, rows, start, size, start + direction * distance)) {
      move(solution, rows, start, size, start - direction * distance);
    }
  }

  private List<Integer> rowsOfBlock(MsaSolution solution, int chosen, int start, int size) {
    boolean wholeBlock = random.nextDouble() < 0.5;
    List<Integer> rows = new ArrayList<>();
    for (int sequence = 0; sequence < solution.variables().size(); sequence++) {
      if (sequence == chosen
          || (hasGaps(solution.variables().get(sequence), start, size)
              && (wholeBlock || random.nextDouble() < 0.5))) {
        rows.add(sequence);
      }
    }
    return rows;
  }

  private static boolean hasGaps(GapIntervals row, int start, int size) {
    boolean gaps = true;
    for (int column = start; column < start + size && gaps; column++) {
      gaps = row.isGapAt(column);
    }
    return gaps;
  }

  /** Moves the block in the rows where the position is valid; tells whether any row moved. */
  private static boolean move(
      MsaSolution solution, List<Integer> rows, int start, int size, int position) {
    boolean moved = false;
    for (int sequence : rows) {
      boolean[] reduced =
          TmBoundaries.removeColumns(solution.variables().get(sequence).toMask(), start, size);
      if (position >= 0 && position <= reduced.length) {
        ProteinSequence protein = solution.proteins().get(sequence);
        int residues = TmBoundaries.residuesBefore(reduced)[position];
        if (!TmBoundaries.splitsHelix(protein, residues)) {
          boolean[] mask = TmBoundaries.insertGaps(reduced, position, size);
          solution.variables().set(sequence, GapIntervals.fromMask(mask));
          moved = true;
        }
      }
    }
    return moved;
  }
}
