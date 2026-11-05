package org.transept.evolver;

import java.util.List;
import org.transept.alignment.MsaSolution;
import org.transept.operator.MergeAdjacentGapRunsMutation;
import org.transept.operator.RandomChoiceMutation;
import org.transept.operator.RandomSource;
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

  /** Crossover names, in the order of the parameter spaces. */
  public static final List<String> CROSSOVERS = List.of(TM_SAFE_SINGLE_POINT);

  /** Mutation names, in the order of the parameter spaces. */
  public static final List<String> MUTATIONS =
      List.of(TM_SAFE_GAP_SHIFT, MERGE_GAP_RUNS, RANDOM_CHOICE);

  private MsaOperators() {}

  /**
   * Builds a crossover.
   *
   * @param name crossover name
   * @param probability crossover probability
   * @param random source of random numbers
   * @return the crossover
   * @throws IllegalArgumentException if the name is unknown
   */
  public static CrossoverOperator<MsaSolution> crossover(
      String name, double probability, RandomSource random) {
    if (!TM_SAFE_SINGLE_POINT.equals(name)) {
      throw new IllegalArgumentException("Unknown crossover: " + name + "; valid: " + CROSSOVERS);
    }
    return new TmSafeSinglePointCrossover(probability, random);
  }

  /**
   * Builds a mutation.
   *
   * @param name mutation name
   * @param probability mutation probability
   * @param structures structures used to recognise residues without structural data
   * @param random source of random numbers
   * @return the mutation
   * @throws IllegalArgumentException if the name is unknown
   */
  public static MutationOperator<MsaSolution> mutation(
      String name, double probability, ProteinStructures structures, RandomSource random) {
    return switch (name) {
      case TM_SAFE_GAP_SHIFT ->
          new TmSafeGapShiftMutation(
              probability, structures, TmSafeGapShiftMutation.DEFAULT_MAX_TRIES, random);
      case MERGE_GAP_RUNS -> new MergeAdjacentGapRunsMutation(probability, random);
      case RANDOM_CHOICE ->
          new RandomChoiceMutation(
              probability,
              List.of(
                  mutation(TM_SAFE_GAP_SHIFT, 1.0, structures, random),
                  mutation(MERGE_GAP_RUNS, 1.0, structures, random)),
              1,
              random);
      default ->
          throw new IllegalArgumentException("Unknown mutation: " + name + "; valid: " + MUTATIONS);
    };
  }
}
