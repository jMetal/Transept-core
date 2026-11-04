package org.tmm2structalign.operator;

import java.util.ArrayList;
import java.util.List;
import org.tmm2structalign.alignment.GapIntervals;
import org.tmm2structalign.alignment.MsaSolution;
import org.tmm2structalign.model.ProteinSequence;
import org.tmm2structalign.structure.ProteinStructures;
import org.uma.jmetal.operator.mutation.MutationOperator;

/**
 * Moves a run of gaps of one sequence to a new position without ever breaking a transmembrane
 * helix: the hard constraint of {@link TmSafeSinglePointCrossover} applied to the mutation.
 *
 * <p>With the mutation probability, a sequence with at least one gap run is chosen at random,
 * and one of its runs is removed. A position {@code k} of the remaining row places the run
 * between residues {@code r - 1} and {@code r}, where {@code r} is the number of residues before
 * {@code k}. The position is <b>valid</b> only if those two residues are not both transmembrane,
 * whatever lies at {@code k} (a residue or an existing gap run). Among valid positions, a
 * candidate is accepted by rejection sampling with the preferences of the original operator:
 *
 * <ul>
 *   <li>inside an existing gap run: {@value #ACCEPT_GAP}
 *   <li>before a residue without structural data: {@value #ACCEPT_UNSTRUCTURED}
 *   <li>before any other residue: {@value #ACCEPT_OTHER}
 * </ul>
 *
 * <p>After {@code maxTries} rejections a valid position is chosen uniformly. If no valid
 * position differs from the original one, the solution is left unchanged.
 *
 * <p>Differences with {@link TmAwareGapShiftMutation}, which is kept unchanged: a helix is never
 * broken (the original only discouraged it, with acceptance 0.05, and not at all in proteins
 * without structure, nor when growing a gap run inside a helix); the boundary is judged by both
 * flanking residues, as in the crossover; any sequence with gaps can be mutated (the original
 * required a run of two or more gaps but then moved any run).
 */
public final class TmSafeGapShiftMutation implements MutationOperator<MsaSolution> {

  /** Acceptance probability when the position is inside an existing gap run. */
  public static final double ACCEPT_GAP = 0.60;

  /** Acceptance probability before a residue without structural data. */
  public static final double ACCEPT_UNSTRUCTURED = 0.95;

  /** Acceptance probability before any other residue. */
  public static final double ACCEPT_OTHER = 0.90;

  /** Number of rejection-sampling attempts, as in the original operator. */
  public static final int DEFAULT_MAX_TRIES = 60;

  private final double probability;
  private final ProteinStructures structures;
  private final int maxTries;
  private final RandomSource random;

  /**
   * Creates the operator with the default number of attempts and the global jMetal generator.
   *
   * @param probability probability of mutating a solution
   * @param structures structures used to recognise residues without structural data
   */
  public TmSafeGapShiftMutation(double probability, ProteinStructures structures) {
    this(probability, structures, DEFAULT_MAX_TRIES, RandomSource.jmetal());
  }

  /**
   * Creates the operator.
   *
   * @param probability probability of mutating a solution
   * @param structures structures used to recognise residues without structural data
   * @param maxTries rejection-sampling attempts before choosing a valid position uniformly
   * @param random source of random numbers
   */
  public TmSafeGapShiftMutation(
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
      List<Integer> eligible = sequencesWithGaps(solution);
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
    boolean[] reduced = TmBoundaries.removeColumns(row.toMask(), start, size);
    ProteinSequence protein = solution.proteins().get(sequence);
    int[] residuesBefore = TmBoundaries.residuesBefore(reduced);
    List<Integer> valid = validPositions(protein, residuesBefore, start);
    GapIntervals result = row;
    if (!valid.isEmpty()) {
      int position = choosePosition(protein, sequence, reduced, residuesBefore, valid);
      result = GapIntervals.fromMask(TmBoundaries.insertGaps(reduced, position, size));
    }
    return result;
  }

  /** Positions of the reduced row whose boundary does not split a helix, except the original. */
  private static List<Integer> validPositions(
      ProteinSequence protein, int[] residuesBefore, int originalStart) {
    List<Integer> valid = new ArrayList<>();
    for (int position = 0; position < residuesBefore.length; position++) {
      if (position != originalStart
          && !TmBoundaries.splitsHelix(protein, residuesBefore[position])) {
        valid.add(position);
      }
    }
    return valid;
  }

  private int choosePosition(
      ProteinSequence protein,
      int sequence,
      boolean[] reduced,
      int[] residuesBefore,
      List<Integer> valid) {
    int accepted = -1;
    for (int attempt = 0; attempt < maxTries && accepted < 0; attempt++) {
      int candidate = valid.get(random.nextInt(0, valid.size() - 1));
      double acceptance = acceptance(protein, sequence, reduced, residuesBefore, candidate);
      if (random.nextDouble() < acceptance) {
        accepted = candidate;
      }
    }
    return accepted >= 0 ? accepted : valid.get(random.nextInt(0, valid.size() - 1));
  }

  private double acceptance(
      ProteinSequence protein,
      int sequence,
      boolean[] reduced,
      int[] residuesBefore,
      int position) {
    int residue = residuesBefore[position];
    double acceptance;
    if (position < reduced.length && reduced[position]) {
      acceptance = ACCEPT_GAP;
    } else if (residue < protein.length() && structures.structureIndex(sequence, residue) < 0) {
      acceptance = ACCEPT_UNSTRUCTURED;
    } else {
      acceptance = ACCEPT_OTHER;
    }
    return acceptance;
  }

  private static List<Integer> sequencesWithGaps(MsaSolution solution) {
    List<Integer> eligible = new ArrayList<>();
    for (int sequence = 0; sequence < solution.variables().size(); sequence++) {
      if (solution.variables().get(sequence).runCount() > 0) {
        eligible.add(sequence);
      }
    }
    return eligible;
  }
}
