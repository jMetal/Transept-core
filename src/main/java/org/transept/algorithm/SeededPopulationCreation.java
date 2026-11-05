package org.transept.algorithm;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import org.transept.alignment.GapColumnCompactor;
import org.transept.alignment.GapIntervals;
import org.transept.alignment.Msa;
import org.transept.alignment.MsaSolution;
import org.transept.operator.RandomSource;
import org.transept.operator.TmSafeSinglePointCrossover;
import org.transept.problem.StructuralMsaProblem;
import org.uma.jmetal.component.catalogue.common.solutionscreation.SolutionsCreation;
import org.uma.jmetal.operator.mutation.MutationOperator;

/**
 * Initial population built from seed alignments (e.g. produced by MAFFT, ClustalW...).
 *
 * <p>Seeds are de-duplicated (identical alignments count once), kept in the given order, and
 * truncated to the population size if there are too many: the first ones in the list are kept.
 * The rest of the population is filled with children of the TM-safe crossover applied to random
 * pairs of members, discarding children equivalent to a member (equal once gap-only columns are
 * removed). If the crossover stops producing new alignments, the remaining members are mutated
 * copies of random members. The result has exactly the requested size (review finding M2).
 */
public final class SeededPopulationCreation implements SolutionsCreation<MsaSolution> {

  /** Crossovers tried per missing member before falling back to mutated copies. */
  public static final int DISTINCT_ATTEMPTS_PER_MEMBER = 10;

  /** Mutations applied to a copy of a member when the crossover produces no new alignment. */
  public static final int PERTURBATION_MUTATIONS = PerturbedSeedCreation.DEFAULT_MUTATIONS;

  private final StructuralMsaProblem problem;
  private final List<Msa> seeds;
  private final int populationSize;
  private final RandomSource random;
  private final MutationOperator<MsaSolution> mutation;

  /**
   * Creates the population builder.
   *
   * @param problem problem whose proteins the seeds align
   * @param seeds seed alignments, in order of preference
   * @param populationSize number of solutions to create
   * @param random source of random numbers for parent selection and crossover
   * @param mutation mutation used for the fallback copies; it should mutate with probability 1
   */
  public SeededPopulationCreation(
      StructuralMsaProblem problem,
      List<Msa> seeds,
      int populationSize,
      RandomSource random,
      MutationOperator<MsaSolution> mutation) {
    if (populationSize < 1) {
      throw new IllegalArgumentException("Population size must be positive: " + populationSize);
    }
    if (seeds == null || seeds.isEmpty()) {
      throw new IllegalArgumentException("At least one seed alignment is required");
    }
    if (mutation == null) {
      throw new IllegalArgumentException("A mutation is required");
    }
    this.problem = problem;
    this.seeds = List.copyOf(seeds);
    this.populationSize = populationSize;
    this.random = random;
    this.mutation = mutation;
  }

  /**
   * Returns the seeds that will be used, without duplicates and truncated to the population size.
   *
   * @return distinct seeds in their original order
   */
  public List<Msa> distinctSeeds() {
    Set<List<GapIntervals>> seen = new LinkedHashSet<>();
    List<Msa> distinct = new ArrayList<>();
    for (Msa seed : seeds) {
      if (seen.add(seed.rows()) && distinct.size() < populationSize) {
        distinct.add(seed);
      }
    }
    return distinct;
  }

  @Override
  public List<MsaSolution> create() {
    List<MsaSolution> population = new ArrayList<>(populationSize);
    for (Msa seed : distinctSeeds()) {
      population.add(problem.createSolution(seed));
    }
    if (population.size() < populationSize && population.size() < 2) {
      throw new IllegalArgumentException(
          "At least two distinct seed alignments are needed to fill a population of "
              + populationSize);
    }
    fillWithDistinctChildren(population, new TmSafeSinglePointCrossover(1.0, random));
    return population;
  }

  /**
   * Reorders a population so that its first {@code seedCount} members (the seeds) are spread
   * over equally spaced positions. MOEA/D assigns the i-th solution to the i-th weight vector, so
   * this spreads the seeds along the front instead of placing them all at one end.
   *
   * @param population population whose first members are the seeds
   * @param seedCount number of seeds at the beginning of the population
   * @return a new list with the seeds at positions {@code round(j (N - 1) / (k - 1))}
   */
  public static List<MsaSolution> spreadSeeds(List<MsaSolution> population, int seedCount) {
    int size = population.size();
    MsaSolution[] spread = new MsaSolution[size];
    for (int seed = 0; seed < seedCount; seed++) {
      int position =
          seedCount == 1
              ? (size - 1) / 2
              : (int) Math.round(seed * (size - 1) / (double) (seedCount - 1));
      spread[position] = population.get(seed);
    }
    int next = seedCount;
    for (int position = 0; position < size; position++) {
      if (spread[position] == null) {
        spread[position] = population.get(next++);
      }
    }
    return new ArrayList<>(List.of(spread));
  }

  private void fillWithDistinctChildren(
      List<MsaSolution> population, TmSafeSinglePointCrossover crossover) {
    Set<List<GapIntervals>> seen = new HashSet<>();
    population.forEach(member -> seen.add(compacted(member)));
    int attempts = DISTINCT_ATTEMPTS_PER_MEMBER * (populationSize - population.size());
    while (population.size() < populationSize && attempts > 0) {
      for (MsaSolution child : crossover.execute(randomParents(population))) {
        if (population.size() < populationSize && seen.add(compacted(child))) {
          population.add(child);
        }
      }
      attempts--;
    }
    while (population.size() < populationSize) {
      MsaSolution child = perturb(population.get(random.nextInt(0, population.size() - 1)).copy());
      if (seen.add(compacted(child)) || attempts-- < -populationSize) {
        population.add(child);
      }
    }
  }

  private List<MsaSolution> randomParents(List<MsaSolution> population) {
    int first = random.nextInt(0, population.size() - 1);
    int second = first;
    while (second == first) {
      second = random.nextInt(0, population.size() - 1);
    }
    return List.of(population.get(first), population.get(second));
  }

  private MsaSolution perturb(MsaSolution solution) {
    MsaSolution perturbed = solution;
    for (int step = 0; step < PERTURBATION_MUTATIONS; step++) {
      perturbed = mutation.execute(perturbed);
    }
    return perturbed;
  }

  private static List<GapIntervals> compacted(MsaSolution solution) {
    return GapColumnCompactor.compact(List.copyOf(solution.variables()));
  }
}
