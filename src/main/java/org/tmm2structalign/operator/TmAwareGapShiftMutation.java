package org.tmm2structalign.operator;

import java.util.ArrayList;
import java.util.List;
import org.tmm2structalign.alignment.GapIntervals;
import org.tmm2structalign.alignment.MsaSolution;
import org.tmm2structalign.model.ProteinSequence;
import org.tmm2structalign.structure.ProteinStructures;
import org.uma.jmetal.operator.mutation.MutationOperator;

/**
 * Moves a block of gaps of one sequence to a new position, preferring loops and avoiding
 * transmembrane helices.
 *
 * <p>With the mutation probability: a sequence that has at least one run of two or more gaps is
 * chosen at random; one of its gap runs (of any length, as in the original operator) is removed;
 * a new position {@code k} is proposed uniformly and accepted by rejection sampling with a
 * probability that depends on what currently occupies {@code k} in that sequence:
 *
 * <ul>
 *   <li>a gap: {@value #ACCEPT_GAP}
 *   <li>a residue without structure: {@value #ACCEPT_UNSTRUCTURED}
 *   <li>a transmembrane residue: {@value #ACCEPT_TRANSMEMBRANE}
 *   <li>any other residue: {@value #ACCEPT_OTHER}
 * </ul>
 *
 * <p>After {@code maxTries} rejections, a uniform position is used. The run is reinserted at
 * {@code k}; the alignment length does not change.
 *
 * <p>Port of {@code BioShiftClosedGapsMSAMutation}. Unlike the original, it does nothing when no
 * sequence has a run of two or more gaps, instead of looping forever (review finding M1).
 */
public final class TmAwareGapShiftMutation implements MutationOperator<MsaSolution> {

  /** Acceptance probability when the proposed position is a gap. */
  public static final double ACCEPT_GAP = 0.60;

  /** Acceptance probability for a residue without structural data. */
  public static final double ACCEPT_UNSTRUCTURED = 0.95;

  /** Acceptance probability for a transmembrane residue. */
  public static final double ACCEPT_TRANSMEMBRANE = 0.05;

  /** Acceptance probability for any other residue. */
  public static final double ACCEPT_OTHER = 0.90;

  /** Number of rejection-sampling attempts used by the original operator. */
  public static final int DEFAULT_MAX_TRIES = 60;

  private final double probability;
  private final ProteinStructures structures;
  private final int maxTries;
  private final RandomSource random;

  /**
   * Creates the operator with the original settings and the global jMetal generator.
   *
   * @param probability probability of mutating a solution
   * @param structures structures used to recognise residues without structural data
   */
  public TmAwareGapShiftMutation(double probability, ProteinStructures structures) {
    this(probability, structures, DEFAULT_MAX_TRIES, RandomSource.jmetal());
  }

  /**
   * Creates the operator.
   *
   * @param probability probability of mutating a solution
   * @param structures structures used to recognise residues without structural data
   * @param maxTries rejection-sampling attempts before falling back to a uniform position
   * @param random source of random numbers
   */
  public TmAwareGapShiftMutation(
      double probability, ProteinStructures structures, int maxTries, RandomSource random) {
    if (!(probability >= 0.0 && probability <= 1.0)) {
      throw new IllegalArgumentException("Mutation probability must be in [0, 1]: " + probability);
    }
    if (maxTries < 1) {
      throw new IllegalArgumentException("maxTries must be positive: " + maxTries);
    }
    this.probability = probability;
    this.structures = structures;
    this.maxTries = maxTries;
    this.random = random;
  }

  @Override
  public MsaSolution execute(MsaSolution solution) {
    if (solution == null) {
      throw new IllegalArgumentException("Solution cannot be null");
    }
    if (random.nextDouble() < probability) {
      List<Integer> eligible = sequencesWithGapBlocks(solution);
      if (!eligible.isEmpty()) {
        int sequence = eligible.get(random.nextInt(0, eligible.size() - 1));
        GapIntervals row = solution.variables().get(sequence);
        int run = random.nextInt(0, row.runCount() - 1);
        solution.variables().set(sequence, shiftRun(solution, sequence, row, run));
      }
    }
    return solution;
  }

  @Override
  public double mutationProbability() {
    return probability;
  }

  private GapIntervals shiftRun(MsaSolution solution, int sequence, GapIntervals row, int run) {
    int start = row.start(run);
    int size = row.runLength(run);
    boolean[] reduced = removeColumns(row.toMask(), start, size);
    int position = proposePosition(solution.proteins().get(sequence), sequence, reduced, start);
    return GapIntervals.fromMask(insertGaps(reduced, position, size));
  }

  private int proposePosition(
      ProteinSequence protein, int sequence, boolean[] reducedGaps, int forbidden) {
    int maxPosition = reducedGaps.length;
    int accepted = -1;
    for (int attempt = 0; attempt < maxTries && accepted < 0; attempt++) {
      int candidate = random.nextInt(0, maxPosition);
      if (candidate != forbidden
          && random.nextDouble() < acceptance(protein, sequence, reducedGaps, candidate)) {
        accepted = candidate;
      }
    }
    while (accepted < 0) {
      int candidate = random.nextInt(0, maxPosition);
      accepted = candidate != forbidden || maxPosition == 0 ? candidate : -1;
    }
    return accepted;
  }

  private double acceptance(
      ProteinSequence protein, int sequence, boolean[] reducedGaps, int column) {
    int residue = residueAt(reducedGaps, column);
    double acceptance;
    if (residue < 0 || residue >= protein.length()) {
      acceptance = ACCEPT_GAP;
    } else if (structures.structureIndex(sequence, residue) < 0) {
      acceptance = ACCEPT_UNSTRUCTURED;
    } else if (protein.isTransmembrane(residue)) {
      acceptance = ACCEPT_TRANSMEMBRANE;
    } else {
      acceptance = ACCEPT_OTHER;
    }
    return acceptance;
  }

  /** Residue index at a column of a gap mask, -1 for a gap or a column past the end. */
  private static int residueAt(boolean[] gaps, int column) {
    int residue = -1;
    if (column < gaps.length && !gaps[column]) {
      residue = 0;
      for (int index = 0; index < column; index++) {
        residue += gaps[index] ? 0 : 1;
      }
    }
    return residue;
  }

  private static List<Integer> sequencesWithGapBlocks(MsaSolution solution) {
    List<Integer> eligible = new ArrayList<>();
    for (int sequence = 0; sequence < solution.variables().size(); sequence++) {
      GapIntervals row = solution.variables().get(sequence);
      for (int run = 0; run < row.runCount(); run++) {
        if (row.runLength(run) > 1) {
          eligible.add(sequence);
          break;
        }
      }
    }
    return eligible;
  }

  private static boolean[] removeColumns(boolean[] mask, int start, int count) {
    boolean[] result = new boolean[mask.length - count];
    System.arraycopy(mask, 0, result, 0, start);
    System.arraycopy(mask, start + count, result, start, mask.length - start - count);
    return result;
  }

  private static boolean[] insertGaps(boolean[] mask, int position, int count) {
    boolean[] result = new boolean[mask.length + count];
    System.arraycopy(mask, 0, result, 0, position);
    for (int index = position; index < position + count; index++) {
      result[index] = true;
    }
    System.arraycopy(mask, position, result, position + count, mask.length - position);
    return result;
  }
}
