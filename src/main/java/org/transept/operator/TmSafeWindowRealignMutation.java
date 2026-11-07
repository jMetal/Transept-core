package org.transept.operator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.transept.alignment.GapIntervals;
import org.transept.alignment.MsaSolution;
import org.transept.model.ProteinSequence;
import org.transept.objective.GapPenalties;
import org.transept.objective.SubstitutionMatrix;
import org.transept.objective.TopologyAwareSumOfPairs;
import org.uma.jmetal.operator.mutation.MutationOperator;

/**
 * Realigns one row inside a window of columns against the other rows, by dynamic programming,
 * without ever breaking a transmembrane helix.
 *
 * <p>With the mutation probability, a row with gaps and one of its gap runs are chosen at random,
 * and then a window of columns that contains the run. The residues of the row in the window are
 * placed again in the columns of the window: their order is kept, the other rows do not change and
 * no column is added, so the alignment keeps its length. The placement maximises the sum of pairs
 * of the row with the rest, with the matrices and gap penalties of {@link TopologyAwareSumOfPairs},
 * and gaps cannot separate two transmembrane residues. If there is no such placement, the solution
 * is left unchanged.
 *
 * <p>The objective inside the window is only the sequence: it does not use the structures. It
 * follows the iterative refinement that realigns one sequence with the rest ("remove first",
 * Wallace et al., 2005), and the local search of H4MSA, which realigns a window (Rubio-Largo et
 * al., 2016).
 */
public final class TmSafeWindowRealignMutation implements MutationOperator<MsaSolution> {

  /** Default number of columns of the window. */
  public static final int DEFAULT_WINDOW = 50;

  private static final double IMPOSSIBLE = Double.NEGATIVE_INFINITY;
  private static final char GAP = '-';

  private final double probability;
  private final int window;
  private final RandomSource random;

  /**
   * Creates the operator.
   *
   * @param probability probability of mutating a solution
   * @param window number of columns of the window
   * @param random source of random numbers
   */
  public TmSafeWindowRealignMutation(double probability, int window, RandomSource random) {
    if (!(probability >= 0.0 && probability <= 1.0)) {
      throw new IllegalArgumentException("Mutation probability must be in [0, 1]: " + probability);
    }
    if (window < 2) {
      throw new IllegalArgumentException("The window needs at least two columns: " + window);
    }
    this.probability = probability;
    this.window = window;
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
        int first = windowStart(row, run);
        int columns = Math.min(window, row.alignedLength() - first);
        solution.variables().set(sequence, realign(solution, sequence, first, columns));
      }
    }
    return solution;
  }

  @Override
  public double mutationProbability() {
    return probability;
  }

  /** First column of a window that contains the run or, if it is longer, its beginning. */
  private int windowStart(GapIntervals row, int run) {
    int earliest = Math.max(0, row.end(run) - window + 1);
    int latest = Math.max(0, Math.min(row.start(run), row.alignedLength() - window));
    return random.nextInt(Math.min(earliest, latest), latest);
  }

  /**
   * Places the residues of a row again inside a window.
   *
   * @param solution the alignment
   * @param sequence the row
   * @param first first column of the window
   * @param columns number of columns of the window
   * @return the new row, or the same one if the residues cannot be placed
   */
  static GapIntervals realign(MsaSolution solution, int sequence, int first, int columns) {
    GapIntervals row = solution.variables().get(sequence);
    ProteinSequence protein = solution.proteins().get(sequence);
    boolean[] mask = row.toMask();
    int base = TmBoundaries.residuesBefore(mask)[first];
    int residues = 0;
    for (int column = first; column < first + columns; column++) {
      residues += mask[column] ? 0 : 1;
    }
    GapIntervals result = row;
    if (residues > 0 && residues < columns) {
      Profile profile = Profile.of(solution, sequence, first, columns);
      boolean[] placement = place(protein, profile, base, residues);
      if (placement != null) {
        boolean[] changed = mask.clone();
        System.arraycopy(placement, 0, changed, first, columns);
        result = GapIntervals.fromMask(changed);
      }
    }
    return result;
  }

  /**
   * Finds the best placement by dynamic programming over two tables, with {@code k} residues placed
   * after {@code c} columns: the last column holds a residue, or a gap.
   *
   * @return gap mask of the window, or {@code null} if no placement respects the helices
   */
  private static boolean[] place(ProteinSequence protein, Profile profile, int base, int residues) {
    int columns = profile.columns();
    double[][] withResidue = filled(residues + 1, columns + 1);
    double[][] withGap = filled(residues + 1, columns + 1);
    withResidue[0][0] = 0.0;
    for (int k = 0; k <= residues; k++) {
      boolean canHaveGap = !TmBoundaries.splitsHelix(protein, base + k);
      for (int c = 1; c <= columns; c++) {
        if (k > 0) {
          withResidue[k][c] =
              Math.max(withResidue[k - 1][c - 1], withGap[k - 1][c - 1])
                  + profile.residueScore(protein, base + k - 1, c - 1);
        }
        if (canHaveGap) {
          double opening = withResidue[k][c - 1] - profile.gapOpening(c - 1);
          withGap[k][c] = Math.max(opening, withGap[k][c - 1]) + profile.gapScore(c - 1);
        }
      }
    }
    return traceBack(profile, residues, withResidue, withGap);
  }

  private static boolean[] traceBack(
      Profile profile, int residues, double[][] withResidue, double[][] withGap) {
    int columns = profile.columns();
    if (Math.max(withResidue[residues][columns], withGap[residues][columns]) == IMPOSSIBLE) {
      return null;
    }
    boolean[] gaps = new boolean[columns];
    int k = residues;
    boolean inGap = withGap[k][columns] > withResidue[k][columns];
    for (int c = columns; c > 0; c--) {
      gaps[c - 1] = inGap;
      if (inGap) {
        inGap = withGap[k][c - 1] > withResidue[k][c - 1] - profile.gapOpening(c - 1);
      } else {
        k--;
        inGap = withGap[k][c - 1] > withResidue[k][c - 1];
      }
    }
    return gaps;
  }

  private static double[][] filled(int rows, int columns) {
    double[][] table = new double[rows][columns];
    for (double[] line : table) {
      Arrays.fill(line, IMPOSSIBLE);
    }
    return table;
  }

  /** What the other rows contribute to each column of the window. */
  private record Profile(
      char[][] residues, boolean[][] transmembrane, GapPenalties tmGaps, GapPenalties otherGaps) {

    private static final SubstitutionMatrix TM_MATRIX = SubstitutionMatrix.PHAT;
    private static final SubstitutionMatrix OTHER_MATRIX = SubstitutionMatrix.BLOSUM62;

    static Profile of(MsaSolution solution, int sequence, int first, int columns) {
      int others = solution.variables().size() - 1;
      char[][] residues = new char[others][columns];
      boolean[][] transmembrane = new boolean[others][columns];
      int line = 0;
      for (int other = 0; other <= others; other++) {
        if (other != sequence) {
          ProteinSequence protein = solution.proteins().get(other);
          int[] residueAt = solution.variables().get(other).columnToResidue();
          for (int c = 0; c < columns; c++) {
            int residue = residueAt[first + c];
            residues[line][c] = residue < 0 ? GAP : protein.residueAt(residue);
            transmembrane[line][c] = residue >= 0 && protein.isTransmembrane(residue);
          }
          line++;
        }
      }
      return new Profile(
          residues,
          transmembrane,
          TopologyAwareSumOfPairs.TRANSMEMBRANE_GAPS,
          TopologyAwareSumOfPairs.GENERAL_GAPS);
    }

    int columns() {
      return residues.length == 0 ? 0 : residues[0].length;
    }

    double residueScore(ProteinSequence protein, int residue, int column) {
      char mine = protein.residueAt(residue);
      boolean mineTm = protein.isTransmembrane(residue);
      double score = 0.0;
      for (int line = 0; line < residues.length; line++) {
        char other = residues[line][column];
        if (other == GAP) {
          score -= (mineTm ? tmGaps : otherGaps).extend();
        } else {
          boolean tm = mineTm || transmembrane[line][column];
          score += (tm ? TM_MATRIX : OTHER_MATRIX).score(mine, other);
        }
      }
      return score;
    }

    double gapScore(int column) {
      double score = 0.0;
      for (int line = 0; line < residues.length; line++) {
        if (residues[line][column] != GAP) {
          score -= (transmembrane[line][column] ? tmGaps : otherGaps).extend();
        }
      }
      return score;
    }

    double gapOpening(int column) {
      double penalty = 0.0;
      for (int line = 0; line < residues.length; line++) {
        if (residues[line][column] != GAP) {
          GapPenalties gaps = transmembrane[line][column] ? tmGaps : otherGaps;
          penalty += gaps.open() - gaps.extend();
        }
      }
      return penalty;
    }
  }
}
