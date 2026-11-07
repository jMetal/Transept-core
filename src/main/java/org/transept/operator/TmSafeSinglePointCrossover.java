package org.transept.operator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import org.transept.alignment.GapIntervals;
import org.transept.alignment.MsaSolution;
import org.transept.model.ProteinSequence;
import org.transept.model.ProteinSet;
import org.uma.jmetal.operator.crossover.CrossoverOperator;

/**
 * Single-point crossover that never cuts a transmembrane helix.
 *
 * <p>A cut is an aligned column {@code c} of the first parent. For every sequence it leaves
 * {@code r_i} residues on the left; the cut is valid only if no sequence has two consecutive
 * transmembrane residues at positions {@code r_i - 1} and {@code r_i}. The child takes columns
 * {@code [0, c]} of the first parent and, for every sequence, the rest of its residues as aligned
 * in the second parent. Sequences whose remaining residues start later in the second parent are
 * padded with gaps so that the second parent's columns stay aligned. If no valid cut exists the
 * child is a copy of the first parent.
 *
 * <p>Port of {@code BioSPXMSACrossover} (without its unused {@code alpha} parameter). Every
 * crossover produces two children, (first | second) and (second | first), each with its own
 * probability draw.
 *
 * <p>The padding depends on the cut: it can be zero, or displace whole blocks of the second parent
 * in alignments with long variable regions. With {@code cutCandidates} above 1, that many valid
 * cuts are drawn and the one that needs the least padding is used.
 */
public final class TmSafeSinglePointCrossover implements CrossoverOperator<MsaSolution> {

  /** Number of candidate cuts that leaves the choice of the cut at random. */
  public static final int DEFAULT_CUT_CANDIDATES = 1;

  private final double probability;
  private final int cutCandidates;
  private final RandomSource random;

  /**
   * Creates the operator using the global jMetal generator.
   *
   * @param probability probability of applying the crossover to each child
   */
  public TmSafeSinglePointCrossover(double probability) {
    this(probability, RandomSource.jmetal());
  }

  /**
   * Creates the operator.
   *
   * @param probability probability of applying the crossover to each child
   * @param random source of random numbers
   */
  public TmSafeSinglePointCrossover(double probability, RandomSource random) {
    this(probability, DEFAULT_CUT_CANDIDATES, random);
  }

  /**
   * Creates the operator.
   *
   * @param probability probability of applying the crossover to each child
   * @param cutCandidates valid cuts drawn per crossover, of which the one with the least padding is
   *     used
   * @param random source of random numbers
   */
  public TmSafeSinglePointCrossover(double probability, int cutCandidates, RandomSource random) {
    if (!(probability >= 0.0 && probability <= 1.0)) {
      throw new IllegalArgumentException("Crossover probability must be in [0, 1]: " + probability);
    }
    if (cutCandidates < 1) {
      throw new IllegalArgumentException("cutCandidates must be positive: " + cutCandidates);
    }
    this.probability = probability;
    this.cutCandidates = cutCandidates;
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
    return List.of(offspring(first, second), offspring(second, first));
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
   * Returns every column of the first parent at which a cut is valid.
   *
   * @param first parent providing the prefix
   * @return valid cut columns in increasing order
   */
  public List<Integer> compatibleCuts(MsaSolution first) {
    List<GapIntervals> rows = first.variables();
    int length = first.alignmentLength();
    int[][] residuesUpTo = residuesUpToColumn(rows);
    List<Integer> cuts = new ArrayList<>();
    for (int cut = 1; cut < length - 1; cut++) {
      if (!isGapOnly(rows, cut) && !splitsHelix(first.proteins(), residuesUpTo, cut)) {
        cuts.add(cut);
      }
    }
    return cuts;
  }

  /**
   * Builds the child of a cut.
   *
   * @param first parent providing columns {@code [0, cut]}
   * @param second parent providing the remaining residues of every sequence
   * @param cut cut column of the first parent
   * @return the child, not evaluated
   */
  public MsaSolution recombine(MsaSolution first, MsaSolution second, int cut) {
    int sequences = first.variables().size();
    int[] suffixStart = suffixStarts(first, second, cut);
    int minimumSuffixStart = Arrays.stream(suffixStart).min().orElse(0);
    List<GapIntervals> childRows = new ArrayList<>(sequences);
    for (int sequence = 0; sequence < sequences; sequence++) {
      boolean[] prefix = first.variables().get(sequence).toMask();
      boolean[] suffix = second.variables().get(sequence).toMask();
      int padding = suffixStart[sequence] - minimumSuffixStart;
      boolean[] child = new boolean[cut + 1 + padding + suffix.length - suffixStart[sequence]];
      System.arraycopy(prefix, 0, child, 0, cut + 1);
      Arrays.fill(child, cut + 1, cut + 1 + padding, true);
      System.arraycopy(
          suffix,
          suffixStart[sequence],
          child,
          cut + 1 + padding,
          suffix.length - suffixStart[sequence]);
      childRows.add(GapIntervals.fromMask(child));
    }
    return new MsaSolution(first.proteins(), childRows, first.objectives().length);
  }

  private MsaSolution offspring(MsaSolution first, MsaSolution second) {
    MsaSolution child = first.copy();
    if (random.nextDouble() < probability) {
      Optional<Integer> cut = chooseCut(first, second);
      if (cut.isPresent()) {
        child = recombine(first, second, cut.get());
      }
    }
    return child;
  }

  private Optional<Integer> chooseCut(MsaSolution first, MsaSolution second) {
    List<Integer> cuts = compatibleCuts(first);
    Optional<Integer> best = Optional.empty();
    long leastPadding = Long.MAX_VALUE;
    for (int draw = 0; draw < cutCandidates && !cuts.isEmpty(); draw++) {
      int cut = cuts.get(random.nextInt(0, cuts.size() - 1));
      long padding = cutCandidates == 1 ? 0 : padding(first, second, cut);
      if (padding < leastPadding) {
        leastPadding = padding;
        best = Optional.of(cut);
      }
    }
    return best;
  }

  /**
   * Returns the gap columns that a cut adds to the child, over all its rows.
   *
   * @param first parent providing columns {@code [0, cut]}
   * @param second parent providing the remaining residues of every sequence
   * @param cut cut column of the first parent
   * @return the total padding
   */
  public static long padding(MsaSolution first, MsaSolution second, int cut) {
    int[] suffixStart = suffixStarts(first, second, cut);
    int minimum = Arrays.stream(suffixStart).min().orElse(0);
    return Arrays.stream(suffixStart).asLongStream().map(start -> start - minimum).sum();
  }

  /** First column of the second parent with the residues that each row takes from it. */
  private static int[] suffixStarts(MsaSolution first, MsaSolution second, int cut) {
    int[] suffixStart = new int[first.variables().size()];
    for (int sequence = 0; sequence < suffixStart.length; sequence++) {
      int residuesOnLeft = residuesUpTo(first.variables().get(sequence), cut);
      suffixStart[sequence] = columnAfterResidues(second.variables().get(sequence), residuesOnLeft);
    }
    return suffixStart;
  }

  private static boolean isGapOnly(List<GapIntervals> rows, int column) {
    return rows.stream().allMatch(row -> row.isGapAt(column));
  }

  private static boolean splitsHelix(ProteinSet proteins, int[][] residuesUpTo, int cut) {
    boolean splits = false;
    for (int sequence = 0; sequence < residuesUpTo.length && !splits; sequence++) {
      ProteinSequence protein = proteins.get(sequence);
      int left = residuesUpTo[sequence][cut];
      splits =
          left > 0
              && left < protein.length()
              && protein.isTransmembrane(left - 1)
              && protein.isTransmembrane(left);
    }
    return splits;
  }

  /** Number of residues of every row in columns {@code [0, column]}. */
  private static int[][] residuesUpToColumn(List<GapIntervals> rows) {
    int[][] counts = new int[rows.size()][];
    for (int sequence = 0; sequence < rows.size(); sequence++) {
      boolean[] gaps = rows.get(sequence).toMask();
      counts[sequence] = new int[gaps.length];
      int residues = 0;
      for (int column = 0; column < gaps.length; column++) {
        residues += gaps[column] ? 0 : 1;
        counts[sequence][column] = residues;
      }
    }
    return counts;
  }

  private static int residuesUpTo(GapIntervals row, int column) {
    int residues = 0;
    for (int index = 0; index <= column; index++) {
      residues += row.isGapAt(index) ? 0 : 1;
    }
    return residues;
  }

  /** First column after the given number of residues of a row (0 if no residue is consumed). */
  private static int columnAfterResidues(GapIntervals row, int residues) {
    int column = 0;
    if (residues > 0) {
      int[] residueAtColumn = row.columnToResidue();
      while (residueAtColumn[column] != residues - 1) {
        column++;
      }
      column++;
    }
    return column;
  }
}
