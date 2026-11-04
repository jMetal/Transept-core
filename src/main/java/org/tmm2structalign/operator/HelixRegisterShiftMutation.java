package org.tmm2structalign.operator;

import java.util.ArrayList;
import java.util.List;
import org.tmm2structalign.alignment.GapIntervals;
import org.tmm2structalign.alignment.MsaSolution;
import org.tmm2structalign.model.ProteinSequence;
import org.uma.jmetal.operator.mutation.MutationOperator;

/**
 * Moves a whole transmembrane helix of one row a few columns to the left or to the right, to
 * correct register errors (a helix aligned one or two positions out of phase with the others).
 *
 * <p>A helix is a maximal run of transmembrane residues. With the mutation probability, one
 * (row, helix) pair is chosen at random among the helices whose residues occupy consecutive
 * columns, together with a direction and a shift {@code s} in {@code [1, maxShift]}. To move the
 * helix left, the {@code s} gap columns of the loop on its left that are closest to it are
 * removed and {@code s} gaps are inserted right after the helix; moving it right is symmetric.
 * The loop is the stretch between the helix and the previous (or next) helix of the row, or the
 * end of the row. The helix and the loop residues between it and the removed gaps shift by
 * {@code s} columns; the rest of the alignment and its length do not change.
 *
 * <p>The helix is never broken: its residues stay together, and the inserted gaps sit between
 * its first or last residue and a loop residue, a valid boundary because the helix is maximal.
 * A direction is only possible when the loop holds at least one gap; if no helix of any row can
 * move, the solution is left unchanged.
 */
public final class HelixRegisterShiftMutation implements MutationOperator<MsaSolution> {

  /** Default maximum number of columns a helix is moved. */
  public static final int DEFAULT_MAX_SHIFT = 3;

  private final double probability;
  private final int maxShift;
  private final RandomSource random;

  /**
   * Creates the operator with the default maximum shift and the global jMetal generator.
   *
   * @param probability probability of mutating a solution
   */
  public HelixRegisterShiftMutation(double probability) {
    this(probability, DEFAULT_MAX_SHIFT, RandomSource.jmetal());
  }

  /**
   * Creates the operator.
   *
   * @param probability probability of mutating a solution
   * @param maxShift maximum number of columns a helix is moved
   * @param random source of random numbers
   */
  public HelixRegisterShiftMutation(double probability, int maxShift, RandomSource random) {
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
      List<Candidate> candidates = new ArrayList<>();
      for (int sequence = 0; sequence < solution.variables().size(); sequence++) {
        candidates.addAll(
            candidates(sequence, solution.proteins().get(sequence),
                solution.variables().get(sequence)));
      }
      if (!candidates.isEmpty()) {
        Candidate chosen = candidates.get(random.nextInt(0, candidates.size() - 1));
        boolean left =
            !chosen.leftGaps().isEmpty()
                && (chosen.rightGaps().isEmpty() || random.nextDouble() < 0.5);
        List<Integer> gaps = left ? chosen.leftGaps() : chosen.rightGaps();
        int shift = random.nextInt(1, Math.min(maxShift, gaps.size()));
        GapIntervals row = solution.variables().get(chosen.sequence());
        solution.variables().set(chosen.sequence(), shift(row, chosen, left, shift));
      }
    }
    return solution;
  }

  @Override
  public double mutationProbability() {
    return probability;
  }

  /** Returns the maximum number of columns a helix is moved. */
  public int maxShift() {
    return maxShift;
  }

  /**
   * A helix of one row that can move.
   *
   * @param sequence row index
   * @param firstColumn column of the first residue of the helix
   * @param lastColumn column of the last residue of the helix
   * @param leftGaps gap columns of the loop on the left, closest to the helix first
   * @param rightGaps gap columns of the loop on the right, closest to the helix first
   */
  record Candidate(
      int sequence,
      int firstColumn,
      int lastColumn,
      List<Integer> leftGaps,
      List<Integer> rightGaps) {}

  /**
   * Finds the helices of a row that occupy consecutive columns and have a gap in at least one of
   * their flanking loops.
   *
   * @param sequence row index
   * @param protein protein of the row
   * @param row gap layout of the row
   * @return the helices that can move, from left to right
   */
  static List<Candidate> candidates(int sequence, ProteinSequence protein, GapIntervals row) {
    boolean[] gaps = row.toMask();
    int[] columnOf = residueColumns(row, protein.length());
    List<int[]> helices = helices(protein);
    List<Candidate> candidates = new ArrayList<>();
    for (int helix = 0; helix < helices.size(); helix++) {
      int first = helices.get(helix)[0];
      int last = helices.get(helix)[1];
      if (columnOf[last] - columnOf[first] == last - first) {
        int leftLimit = helix > 0 ? columnOf[helices.get(helix - 1)[1]] + 1 : 0;
        int rightLimit =
            helix + 1 < helices.size() ? columnOf[helices.get(helix + 1)[0]] - 1 : gaps.length - 1;
        List<Integer> leftGaps = new ArrayList<>();
        for (int column = columnOf[first] - 1; column >= leftLimit; column--) {
          if (gaps[column]) {
            leftGaps.add(column);
          }
        }
        List<Integer> rightGaps = new ArrayList<>();
        for (int column = columnOf[last] + 1; column <= rightLimit; column++) {
          if (gaps[column]) {
            rightGaps.add(column);
          }
        }
        if (!leftGaps.isEmpty() || !rightGaps.isEmpty()) {
          candidates.add(
              new Candidate(
                  sequence,
                  columnOf[first],
                  columnOf[last],
                  List.copyOf(leftGaps),
                  List.copyOf(rightGaps)));
        }
      }
    }
    return candidates;
  }

  /**
   * Moves a helix.
   *
   * @param row gap layout of the row
   * @param helix the helix to move
   * @param left true to move it left, false to move it right
   * @param shift number of columns, at most the number of gaps of the loop in that direction
   * @return the new gap layout
   */
  static GapIntervals shift(GapIntervals row, Candidate helix, boolean left, int shift) {
    boolean[] gaps = row.toMask();
    List<Integer> removed = (left ? helix.leftGaps() : helix.rightGaps()).subList(0, shift);
    List<Integer> descending = new ArrayList<>(removed);
    descending.sort((first, second) -> Integer.compare(second, first));
    for (int column : descending) {
      gaps = TmBoundaries.removeColumns(gaps, column, 1);
    }
    int insertion = left ? helix.lastColumn() - shift + 1 : helix.firstColumn();
    return GapIntervals.fromMask(TmBoundaries.insertGaps(gaps, insertion, shift));
  }

  /** Maximal runs of transmembrane residues, as inclusive {@code [first, last]} pairs. */
  private static List<int[]> helices(ProteinSequence protein) {
    List<int[]> helices = new ArrayList<>();
    int residue = 0;
    while (residue < protein.length()) {
      if (protein.isTransmembrane(residue)) {
        int first = residue;
        while (residue + 1 < protein.length() && protein.isTransmembrane(residue + 1)) {
          residue++;
        }
        helices.add(new int[] {first, residue});
      }
      residue++;
    }
    return helices;
  }

  private static int[] residueColumns(GapIntervals row, int residues) {
    int[] columnOf = new int[residues];
    int[] residueAtColumn = row.columnToResidue();
    for (int column = 0; column < residueAtColumn.length; column++) {
      if (residueAtColumn[column] >= 0) {
        columnOf[residueAtColumn[column]] = column;
      }
    }
    return columnOf;
  }
}
