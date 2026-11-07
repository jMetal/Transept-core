package org.transept.algorithm;

import java.time.Duration;
import org.transept.operator.TmSafeGapShiftMutation;

/**
 * Parameters of the evolutionary algorithm.
 *
 * @param maxEvaluations stopping criterion, number of evaluations (ignored if {@code maxTime}
 *     is set)
 * @param populationSize population size
 * @param offspringPopulationSize offspring created per generation
 * @param crossoverProbability probability of the TM-safe crossover
 * @param mutationProbability probability of the gap-shift mutation
 * @param threads number of evaluation threads; 1 evaluates sequentially
 * @param mutationType mutation operator
 * @param maxTime stopping criterion by computing time, or {@code null} to stop by evaluations
 * @param configuration configuration string ({@code --parameter value ...}) of the algorithms
 *     built from a parameter space; empty to use their default configuration
 * @param maxShift largest number of columns the mutation can move a gap run; 0 for no limit
 */
public record AlgorithmConfiguration(
    int maxEvaluations,
    int populationSize,
    int offspringPopulationSize,
    double crossoverProbability,
    double mutationProbability,
    int threads,
    MutationType mutationType,
    Duration maxTime,
    String configuration,
    int maxShift) {

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
    if (mutationType == null) {
      throw new IllegalArgumentException("mutationType is required");
    }
    if (maxTime != null && (maxTime.isZero() || maxTime.isNegative())) {
      throw new IllegalArgumentException("maxTime must be positive: " + maxTime);
    }
    configuration = configuration == null ? "" : configuration.strip();
    if (maxShift < 0) {
      throw new IllegalArgumentException("maxShift cannot be negative: " + maxShift);
    }
  }

  /**
   * Creates a configuration with the default settings of the operators.
   *
   * @param maxEvaluations stopping criterion, number of evaluations (ignored if {@code maxTime}
   *     is set)
   * @param populationSize population size
   * @param offspringPopulationSize offspring created per generation
   * @param crossoverProbability probability of the TM-safe crossover
   * @param mutationProbability probability of the gap-shift mutation
   * @param threads number of evaluation threads; 1 evaluates sequentially
   * @param mutationType mutation operator
   * @param maxTime stopping criterion by computing time, or {@code null} to stop by evaluations
   * @param configuration configuration string of the algorithms built from a parameter space
   */
  public AlgorithmConfiguration(
      int maxEvaluations,
      int populationSize,
      int offspringPopulationSize,
      double crossoverProbability,
      double mutationProbability,
      int threads,
      MutationType mutationType,
      Duration maxTime,
      String configuration) {
    this(
        maxEvaluations,
        populationSize,
        offspringPopulationSize,
        crossoverProbability,
        mutationProbability,
        threads,
        mutationType,
        maxTime,
        configuration,
        TmSafeGapShiftMutation.UNLIMITED_SHIFT);
  }

  /**
   * Creates a configuration that stops by evaluations, with the default configuration string.
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
        null,
        "");
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
