package org.tmm2structalign.algorithm;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import org.tmm2structalign.alignment.GapIntervals;
import org.tmm2structalign.alignment.Msa;
import org.tmm2structalign.alignment.MsaSolution;
import org.tmm2structalign.operator.RandomSource;
import org.tmm2structalign.operator.TmSafeSinglePointCrossover;
import org.tmm2structalign.problem.StructuralMsaProblem;
import org.uma.jmetal.component.catalogue.common.solutionscreation.SolutionsCreation;

/**
 * Initial population built from seed alignments (e.g. produced by MAFFT, ClustalW...).
 *
 * <p>Seeds are de-duplicated (identical alignments count once), kept in the given order, and
 * truncated to the population size if there are too many. The rest of the population is filled
 * with children of the TM-safe crossover applied to random pairs of distinct members, as in the
 * original implementation. The result has exactly the requested size (review finding M2).
 */
public final class SeededPopulationCreation implements SolutionsCreation<MsaSolution> {

  private final StructuralMsaProblem problem;
  private final List<Msa> seeds;
  private final int populationSize;
  private final RandomSource random;

  /**
   * Creates the population builder.
   *
   * @param problem problem whose proteins the seeds align
   * @param seeds seed alignments, in order of preference
   * @param populationSize number of solutions to create
   * @param random source of random numbers for parent selection and crossover
   */
  public SeededPopulationCreation(
      StructuralMsaProblem problem, List<Msa> seeds, int populationSize, RandomSource random) {
    if (populationSize < 1) {
      throw new IllegalArgumentException("Population size must be positive: " + populationSize);
    }
    if (seeds == null || seeds.isEmpty()) {
      throw new IllegalArgumentException("At least one seed alignment is required");
    }
    this.problem = problem;
    this.seeds = List.copyOf(seeds);
    this.populationSize = populationSize;
    this.random = random;
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
    TmSafeSinglePointCrossover crossover = new TmSafeSinglePointCrossover(1.0, random);
    while (population.size() < populationSize) {
      int first = random.nextInt(0, population.size() - 1);
      int second = first;
      while (second == first) {
        second = random.nextInt(0, population.size() - 1);
      }
      List<MsaSolution> parents = List.of(population.get(first), population.get(second));
      for (MsaSolution child : crossover.execute(parents)) {
        if (population.size() < populationSize) {
          population.add(child);
        }
      }
    }
    return population;
  }
}
