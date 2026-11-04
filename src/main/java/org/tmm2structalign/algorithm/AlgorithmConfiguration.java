package org.tmm2structalign.algorithm;

/**
 * Parameters of the evolutionary algorithm.
 *
 * @param maxEvaluations stopping criterion, number of evaluations
 * @param populationSize population size
 * @param offspringPopulationSize offspring created per generation
 * @param crossoverProbability probability of the TM-safe crossover
 * @param mutationProbability probability of the gap-shift mutation
 * @param threads number of evaluation threads; 1 evaluates sequentially
 * @param mutationType mutation operator
 * @param initializationType how the non-seed members of the initial population are created
 */
public record AlgorithmConfiguration(
    int maxEvaluations,
    int populationSize,
    int offspringPopulationSize,
    double crossoverProbability,
    double mutationProbability,
    int threads,
    MutationType mutationType,
    InitializationType initializationType) {

  /** Crossover probability used in the published experiments. */
  public static final double DEFAULT_CROSSOVER_PROBABILITY = 0.8;

  /** Mutation probability used in the published experiments. */
  public static final double DEFAULT_MUTATION_PROBABILITY = 0.2;

  /** Validates the parameters. */
  public AlgorithmConfiguration {
    requirePositive(maxEvaluations, "maxEvaluations");
    requirePositive(populationSize, "populationSize");
    requirePositive(offspringPopulationSize, "offspringPopulationSize");
    requirePositive(threads, "threads");
    requireProbability(crossoverProbability, "crossoverProbability");
    requireProbability(mutationProbability, "mutationProbability");
    if (mutationType == null || initializationType == null) {
      throw new IllegalArgumentException("mutationType and initializationType are required");
    }
  }

  /**
   * Creates a configuration with the original initialization of the population.
   *
   * @param maxEvaluations stopping criterion, number of evaluations
   * @param populationSize population size
   * @param offspringPopulationSize offspring created per generation
   * @param crossoverProbability probability of the TM-safe crossover
   * @param mutationProbability probability of the gap-shift mutation
   * @param threads number of evaluation threads; 1 evaluates sequentially
   * @param mutationType mutation operator
   */
  public AlgorithmConfiguration(
      int maxEvaluations,
      int populationSize,
      int offspringPopulationSize,
      double crossoverProbability,
      double mutationProbability,
      int threads,
      MutationType mutationType) {
    this(
        maxEvaluations,
        populationSize,
        offspringPopulationSize,
        crossoverProbability,
        mutationProbability,
        threads,
        mutationType,
        InitializationType.CROSSOVER);
  }

  /**
   * Creates a configuration with the default mutation, the TM-safe gap shift.
   *
   * @param maxEvaluations stopping criterion, number of evaluations
   * @param populationSize population size
   * @param offspringPopulationSize offspring created per generation
   * @param crossoverProbability probability of the TM-safe crossover
   * @param mutationProbability probability of the gap-shift mutation
   * @param threads number of evaluation threads; 1 evaluates sequentially
   */
  public AlgorithmConfiguration(
      int maxEvaluations,
      int populationSize,
      int offspringPopulationSize,
      double crossoverProbability,
      double mutationProbability,
      int threads) {
    this(
        maxEvaluations,
        populationSize,
        offspringPopulationSize,
        crossoverProbability,
        mutationProbability,
        threads,
        MutationType.TM_SAFE_GAP_SHIFT);
  }

  /**
   * Creates a configuration with the published operator probabilities and an offspring
   * population equal to the population.
   *
   * @param maxEvaluations stopping criterion
   * @param populationSize population size
   * @param threads number of evaluation threads
   * @return the configuration
   */
  public static AlgorithmConfiguration withDefaults(
      int maxEvaluations, int populationSize, int threads) {
    return new AlgorithmConfiguration(
        maxEvaluations,
        populationSize,
        populationSize,
        DEFAULT_CROSSOVER_PROBABILITY,
        DEFAULT_MUTATION_PROBABILITY,
        threads);
  }

  private static void requirePositive(int value, String name) {
    if (value <= 0) {
      throw new IllegalArgumentException(name + " must be positive: " + value);
    }
  }

  private static void requireProbability(double value, String name) {
    if (!(value >= 0.0 && value <= 1.0)) {
      throw new IllegalArgumentException(name + " must be in [0, 1]: " + value);
    }
  }
}
