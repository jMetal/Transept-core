package org.transept.operator;

import java.util.List;
import org.transept.alignment.MsaSolution;
import org.uma.jmetal.operator.mutation.MutationOperator;

/**
 * Applies mutations chosen uniformly at random from a list.
 *
 * <p>With the mutation probability, {@code mutationsPerApplication} operators are drawn with
 * replacement from the list and applied in sequence. The operators of the list should mutate
 * with probability 1. Uniform choice is as effective as adaptive operator scheduling in SAGA
 * (Thomsen and Boomsma, 2004). It replaces {@code MultipleMSAMutation} of M2Align, which always
 * applied the whole list.
 */
public final class RandomChoiceMutation implements MutationOperator<MsaSolution> {

  private final double probability;
  private final List<MutationOperator<MsaSolution>> operators;
  private final int mutationsPerApplication;
  private final RandomSource random;

  /**
   * Creates the operator.
   *
   * @param probability probability of mutating a solution
   * @param operators mutations to choose from, not empty
   * @param mutationsPerApplication operators applied each time the solution is mutated
   * @param random source of random numbers
   */
  public RandomChoiceMutation(
      double probability,
      List<MutationOperator<MsaSolution>> operators,
      int mutationsPerApplication,
      RandomSource random) {
    if (!(probability >= 0.0 && probability <= 1.0)) {
      throw new IllegalArgumentException("Mutation probability must be in [0, 1]: " + probability);
    }
    if (operators == null || operators.isEmpty()) {
      throw new IllegalArgumentException("At least one mutation operator is required");
    }
    if (mutationsPerApplication < 1) {
      throw new IllegalArgumentException(
          "mutationsPerApplication must be positive: " + mutationsPerApplication);
    }
    this.probability = probability;
    this.operators = List.copyOf(operators);
    this.mutationsPerApplication = mutationsPerApplication;
    this.random = random;
  }

  @Override
  public MsaSolution execute(MsaSolution solution) {
    if (solution == null) {
      throw new IllegalArgumentException("Solution cannot be null");
    }
    MsaSolution result = solution;
    if (random.nextDouble() < probability) {
      for (int mutation = 0; mutation < mutationsPerApplication; mutation++) {
        result = operators.get(random.nextInt(0, operators.size() - 1)).execute(result);
      }
    }
    return result;
  }

  @Override
  public double mutationProbability() {
    return probability;
  }

  /** Returns the mutations to choose from. */
  public List<MutationOperator<MsaSolution>> operators() {
    return operators;
  }
}
