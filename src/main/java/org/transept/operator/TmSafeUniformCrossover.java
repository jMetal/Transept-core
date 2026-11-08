package org.transept.operator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.transept.alignment.GapIntervals;
import org.transept.alignment.MsaSolution;
import org.uma.jmetal.operator.crossover.CrossoverOperator;

/**
 * Uniform crossover between the positions where both parents agree, which never adds a gap.
 *
 * <p>A consistent position is a place where both parents have aligned the same residues on its left
 * in every row: after some column of the first parent and some column of the second one, every
 * row has used the same number of residues. Between two consecutive consistent positions, a block
 * of the first parent and a block of the second hold the same residues of every row, aligned in two
 * different ways, so either block can take the place of the other.
 *
 * <p>With the crossover probability, every block is taken at random from one parent or the other,
 * and the second child takes the blocks the first did not. No padding is needed, so the children
 * keep the parents' alignment inside each block. A child separates two transmembrane residues only
 * where one of its parents already did. If the parents agree only at the start and at the end,
 * the children are copies of the parents.
 *
 * <p>It is the uniform crossover of SAGA (Notredame and Higgins, 1996).
 */
public final class TmSafeUniformCrossover implements CrossoverOperator<MsaSolution> {

  private final double probability;
  private final RandomSource random;

  /**
   * Creates the operator.
   *
   * @param probability probability of crossing the two parents
   * @param random source of random numbers
   */
  public TmSafeUniformCrossover(double probability, RandomSource random) {
    if (!(probability >= 0.0 && probability <= 1.0)) {
      throw new IllegalArgumentException("Crossover probability must be in [0, 1]: " + probability);
    }
    this.probability = probability;
    this.random = random;
  }

  @Override
  public List<MsaSolution> execute(List<MsaSolution> parents) {
    if (parents == null || parents.size() != 2) {
      throw new IllegalArgumentException("The crossover needs exactly two parents");
    }
    MsaSolution first = parents.get(0);
    MsaSolution second = parents.get(1);
    if (first.proteins() != second.proteins()) {
      throw new IllegalArgumentException("Parents must align the same protein set");
    }
    List<MsaSolution> children = List.of(first.copy(), second.copy());
    if (random.nextDouble() < probability) {
      List<int[]> positions = consistentPositions(first, second);
      if (positions.size() > 2) {
        boolean[] fromFirst = new boolean[positions.size() - 1];
        for (int block = 0; block < fromFirst.length; block++) {
          fromFirst[block] = random.nextDouble() < 0.5;
        }
        children =
            List.of(
                child(first, second, positions, fromFirst, true),
                child(first, second, positions, fromFirst, false));
      }
    }
    return children;
  }

  @Override
  public double crossoverProbability() {
    return probability;
  }

  @Override
  public int numberOfRequiredParents() {
    return 2;
  }

  @Override
  public int numberOfGeneratedChildren() {
    return 2;
  }

  /**
   * Returns the consistent positions of two parents, as pairs {number of columns of the first
   * parent on the left, number of columns of the second parent on the left}, from the start to
   * the end of both alignments.
   *
   * @param first a parent
   * @param second the other parent
   * @return the positions, in increasing order
   */
  public static List<int[]> consistentPositions(MsaSolution first, MsaSolution second) {
    Map<Residues, Integer> inSecond = new HashMap<>();
    int[][] secondCounts = residuesOnTheLeft(second);
    for (int columns = secondCounts.length - 1; columns >= 0; columns--) {
      inSecond.put(new Residues(secondCounts[columns]), columns);
    }
    List<int[]> positions = new ArrayList<>();
    Residues previous = null;
    int[][] firstCounts = residuesOnTheLeft(first);
    for (int columns = 0; columns < firstCounts.length; columns++) {
      Residues residues = new Residues(firstCounts[columns]);
      Integer match = inSecond.get(residues);
      if (match != null && !residues.equals(previous)) {
        positions.add(new int[] {columns, match});
        previous = residues;
      }
    }
    return positions;
  }

  /** Residues of every row in the first {@code k} columns, for every {@code k} up to the end. */
  private static int[][] residuesOnTheLeft(MsaSolution solution) {
    List<GapIntervals> rows = solution.variables();
    int length = solution.alignmentLength();
    int[][] counts = new int[length + 1][rows.size()];
    for (int sequence = 0; sequence < rows.size(); sequence++) {
      boolean[] gaps = rows.get(sequence).toMask();
      for (int column = 0; column < length; column++) {
        counts[column + 1][sequence] = counts[column][sequence] + (gaps[column] ? 0 : 1);
      }
    }
    return counts;
  }

  private static MsaSolution child(
      MsaSolution first,
      MsaSolution second,
      List<int[]> positions,
      boolean[] fromFirst,
      boolean firstChild) {
    int length = 0;
    for (int block = 0; block < fromFirst.length; block++) {
      int parent = fromFirst[block] == firstChild ? 0 : 1;
      length += positions.get(block + 1)[parent] - positions.get(block)[parent];
    }
    List<GapIntervals> rows = new ArrayList<>();
    for (int sequence = 0; sequence < first.variables().size(); sequence++) {
      boolean[] a = first.variables().get(sequence).toMask();
      boolean[] b = second.variables().get(sequence).toMask();
      boolean[] row = new boolean[length];
      int filled = 0;
      for (int block = 0; block < fromFirst.length; block++) {
        boolean takeFirst = fromFirst[block] == firstChild;
        int parent = takeFirst ? 0 : 1;
        int start = positions.get(block)[parent];
        int size = positions.get(block + 1)[parent] - start;
        System.arraycopy(takeFirst ? a : b, start, row, filled, size);
        filled += size;
      }
      rows.add(GapIntervals.fromMask(row));
    }
    return new MsaSolution(first.proteins(), rows, first.objectives().length);
  }

  /** The number of residues of every row, compared by value. */
  private record Residues(int[] counts) {
    @Override
    public boolean equals(Object other) {
      return other instanceof Residues residues && Arrays.equals(counts, residues.counts);
    }

    @Override
    public int hashCode() {
      return Arrays.hashCode(counts);
    }

    @Override
    public String toString() {
      return Arrays.toString(counts);
    }
  }
}
