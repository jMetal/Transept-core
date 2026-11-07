package org.transept.evolver;

import java.util.List;
import org.transept.alignment.MsaSolution;
import org.transept.operator.MergeAdjacentGapRunsMutation;
import org.transept.operator.RandomChoiceMutation;
import org.transept.operator.RandomSource;
import org.transept.operator.TmSafeBlockShiftMutation;
import org.transept.operator.TmSafeGapShiftMutation;
import org.transept.operator.TmSafeSinglePointCrossover;
import org.transept.structure.ProteinStructures;
import org.uma.jmetal.operator.crossover.CrossoverOperator;
import org.uma.jmetal.operator.mutation.MutationOperator;

/** Builds the variation operators by the names used in the parameter spaces. */
public final class MsaOperators {

  /** Crossover that never cuts a transmembrane helix. */
  public static final String TM_SAFE_SINGLE_POINT = "tmSafeSinglePoint";

  /** Gap shift that never breaks a transmembrane helix. */
  public static final String TM_SAFE_GAP_SHIFT = "tmSafeGapShift";

  /** Joins two consecutive gap runs of a row. */
  public static final String MERGE_GAP_RUNS = "mergeGapRuns";

  /** {@link #TM_SAFE_GAP_SHIFT} or {@link #MERGE_GAP_RUNS}, chosen at random each time. */
  public static final String RANDOM_CHOICE = "randomChoice";

  /** Moves a block of gaps a few columns in several rows at once. */
  public static final String TM_SAFE_BLOCK_SHIFT = "tmSafeBlockShift";

  /** Crossover names, in the order of the parameter spaces. */
  public static final List<String> CROSSOVERS = List.of(TM_SAFE_SINGLE_POINT);

  /** Mutation names, in the order of the parameter spaces. */
  public static final List<String> MUTATIONS =
      List.of(TM_SAFE_GAP_SHIFT, MERGE_GAP_RUNS, RANDOM_CHOICE, TM_SAFE_BLOCK_SHIFT);

  private MsaOperators() {}

  /**
   * Builds a crossover.
   *
   * @param name crossover name
   * @param probability crossover probability
   * @param cutCandidates valid cuts drawn per crossover, of which the one with the least padding is
   *     used
   * @param random source of random numbers
   * @return the crossover
   * @throws IllegalArgumentException if the name is unknown
   */
  public static CrossoverOperator<MsaSolution> crossover(
      String name, double probability, int cutCandidates, RandomSource random) {
    if (!TM_SAFE_SINGLE_POINT.equals(name)) {
      throw new IllegalArgumentException("Unknown crossover: " + name + "; valid: " + CROSSOVERS);
    }
    return new TmSafeSinglePointCrossover(probability, cutCandidates, random);
  }

  /**
   * Builds a mutation.
   *
   * @param name mutation name
   * @param probability mutation probability
   * @param maxShift largest number of columns a gap run can move; {@link
   *     TmSafeGapShiftMutation#UNLIMITED_SHIFT} for no limit
   * @param structures structures used to recognise residues without structural data
   * @param random source of random numbers
   * @return the mutation
   * @throws IllegalArgumentException if the name is unknown
   */
  public static MutationOperator<MsaSolution> mutation(
      String name,
      double probability,
      int maxShift,
      ProteinStructures structures,
      RandomSource random) {
    int blockShift = TmSafeBlockShiftMutation.DEFAULT_MAX_SHIFT;
    return mutation(name, probability, maxShift, blockShift, structures, random);
  }

  /**
   * Builds a mutation.
   *
   * @param name mutation name
   * @param probability mutation probability
   * @param maxShift largest number of columns a gap run can move; {@link
   *     TmSafeGapShiftMutation#UNLIMITED_SHIFT} for no limit
   * @param blockShift largest number of columns a block of gaps can move
   * @param structures structures used to recognise residues without structural data
   * @param random source of random numbers
   * @return the mutation
   * @throws IllegalArgumentException if the name is unknown
   */
  public static MutationOperator<MsaSolution> mutation(
      String name,
      double probability,
      int maxShift,
      int blockShift,
      ProteinStructures structures,
      RandomSource random) {
    return switch (name) {
      case TM_SAFE_GAP_SHIFT ->
          new TmSafeGapShiftMutation(
              probability, structures, TmSafeGapShiftMutation.DEFAULT_MAX_TRIES, maxShift, random);
      case MERGE_GAP_RUNS -> new MergeAdjacentGapRunsMutation(probability, random);
      case RANDOM_CHOICE ->
          new RandomChoiceMutation(
              probability,
              List.of(
                  mutation(TM_SAFE_GAP_SHIFT, 1.0, maxShift, structures, random),
                  mutation(MERGE_GAP_RUNS, 1.0, maxShift, structures, random)),
              1,
              random);
      case TM_SAFE_BLOCK_SHIFT -> new TmSafeBlockShiftMutation(probability, blockShift, random);
      default ->
          throw new IllegalArgumentException("Unknown mutation: " + name + "; valid: " + MUTATIONS);
    };
  }
}
