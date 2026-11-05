package org.transept.algorithm;

import java.util.ArrayList;
import java.util.List;
import org.transept.alignment.Msa;
import org.transept.alignment.MsaSolution;
import org.transept.operator.RandomSource;
import org.transept.problem.StructuralMsaProblem;
import org.uma.jmetal.component.catalogue.common.solutionscreation.SolutionsCreation;
import org.uma.jmetal.operator.mutation.MutationOperator;

/**
 * Creates alignments by applying a fixed number of mutations to a seed chosen at random, as a
 * strong baseline: random search around the seeds, with no selection or recombination.
 */
public final class PerturbedSeedCreation implements SolutionsCreation<MsaSolution> {

  /** Number of mutations applied to each seed by default. */
  public static final int DEFAULT_MUTATIONS = 10;

  private final StructuralMsaProblem problem;
  private final List<Msa> seeds;
  private final int batchSize;
  private final int mutations;
  private final MutationOperator<MsaSolution> mutation;
  private final RandomSource random;

  /**
   * Creates the generator.
   *
   * @param problem problem whose proteins the seeds align
   * @param seeds seed alignments, at least one
   * @param batchSize number of alignments created per call
   * @param mutations number of mutations applied to each copy of a seed, non-negative
   * @param mutation mutation operator, usually with probability 1
   * @param random source of random numbers for choosing seeds
   */
  public PerturbedSeedCreation(
      StructuralMsaProblem problem,
      List<Msa> seeds,
      int batchSize,
      int mutations,
      MutationOperator<MsaSolution> mutation,
      RandomSource random) {
    if (seeds == null || seeds.isEmpty()) {
      throw new IllegalArgumentException("At least one seed alignment is required");
    }
    if (batchSize < 1) {
      throw new IllegalArgumentException("Batch size must be positive: " + batchSize);
    }
    if (mutations < 0) {
      throw new IllegalArgumentException("Number of mutations cannot be negative: " + mutations);
    }
    this.problem = problem;
    this.seeds = List.copyOf(seeds);
    this.batchSize = batchSize;
    this.mutations = mutations;
    this.mutation = mutation;
    this.random = random;
  }

  @Override
  public List<MsaSolution> create() {
    List<MsaSolution> solutions = new ArrayList<>(batchSize);
    for (int index = 0; index < batchSize; index++) {
      Msa seed = seeds.get(random.nextInt(0, seeds.size() - 1));
      MsaSolution solution = problem.createSolution(seed);
      for (int step = 0; step < mutations; step++) {
        solution = mutation.execute(solution);
      }
      solutions.add(solution);
    }
    return solutions;
  }
}
