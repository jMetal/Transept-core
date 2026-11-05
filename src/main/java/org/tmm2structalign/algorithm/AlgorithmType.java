package org.tmm2structalign.algorithm;

import org.tmm2structalign.evolver.MsaGeneticAlgorithm;
import org.tmm2structalign.evolver.MsaNsga2;
import org.tmm2structalign.evolver.MsaSpea2;

/** Search algorithms available to align a dataset. */
public enum AlgorithmType {
  /** NSGA-II with the published configuration (default); configurable. */
  NSGA_II(false, true, false, MsaNsga2.PARAMETER_SPACE),
  /**
   * Asynchronous NSGA-II: every evaluated solution enters the population at once, so all threads
   * stay busy; runs do not repeat exactly.
   */
  ASYNC_NSGA_II(false, true, false, null),
  /** SPEA2-style algorithm (strength ranking, k-nearest-neighbour density); configurable. */
  SPEA2(false, true, false, MsaSpea2.PARAMETER_SPACE),
  /** NSGA-III; the population size must be a multiple of 4. */
  NSGA_III(false, true, false, null),
  /** AGE-MOEA. */
  AGE_MOEA(false, true, false, null),
  /** AGE-MOEA-II. */
  AGE_MOEA_II(false, true, false, null),
  /** SMS-EMOA, steady state: one offspring per iteration. */
  SMS_EMOA(true, true, false, null),
  /** MOEA/D with normalised PBI aggregation, steady state. */
  MOEA_D(true, true, false, null),
  /** PAES, a (1+1) evolution strategy with a crowding-distance archive; mutation only. */
  PAES(true, false, false, null),
  /**
   * Generational genetic algorithm for a single objective, with the seeded population and the
   * operators of NSGA-II; configurable.
   */
  GENETIC_ALGORITHM(false, true, true, MsaGeneticAlgorithm.PARAMETER_SPACE),
  /** Asynchronous steady-state genetic algorithm for a single objective; runs do not repeat. */
  ASYNC_GENETIC_ALGORITHM(false, true, true, null),
  /** Random search over alignments with gaps placed at random (weak baseline). */
  RANDOM_SEARCH(false, false, false, null),
  /** Random search over seeds perturbed with random mutations (strong baseline). */
  RANDOM_SEEDED(false, false, false, null);

  /** Largest number of objectives of the formulations of the tool. */
  public static final int MAX_OBJECTIVES = 3;

  private final boolean steadyState;
  private final boolean usesCrossover;
  private final boolean singleObjective;
  private final String parameterSpace;

  AlgorithmType(
      boolean steadyState, boolean usesCrossover, boolean singleObjective, String parameterSpace) {
    this.steadyState = steadyState;
    this.usesCrossover = usesCrossover;
    this.singleObjective = singleObjective;
    this.parameterSpace = parameterSpace;
  }

  /** Returns true if one solution is evaluated per iteration, so threads do not help. */
  public boolean isSteadyState() {
    return steadyState;
  }

  /** Returns true if the algorithm applies the crossover operator. */
  public boolean usesCrossover() {
    return usesCrossover;
  }

  /** Returns true if the algorithm optimises a single objective. */
  public boolean isSingleObjective() {
    return singleObjective;
  }

  /** Returns true if the algorithm evaluates asynchronously, so its runs do not repeat exactly. */
  public boolean isAsynchronous() {
    return name().startsWith("ASYNC_");
  }

  /** Returns the smallest number of objectives the algorithm accepts. */
  public int minObjectives() {
    return singleObjective || isRandomSearch() ? 1 : 2;
  }

  /** Returns the largest number of objectives the algorithm accepts. */
  public int maxObjectives() {
    return singleObjective ? 1 : MAX_OBJECTIVES;
  }

  /** Returns true for the random searches, which accept any number of objectives. */
  public boolean isRandomSearch() {
    return this == RANDOM_SEARCH || this == RANDOM_SEEDED;
  }

  /** Returns true if the algorithm is built from a parameter space and accepts a configuration. */
  public boolean isConfigurable() {
    return parameterSpace != null;
  }

  /** Returns the parameter space resource, or {@code null} if the algorithm is not configurable. */
  public String parameterSpace() {
    return parameterSpace;
  }
}
