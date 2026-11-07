package org.transept.evolver;

import java.util.List;
import org.transept.alignment.MsaSolution;
import org.transept.operator.RandomSource;
import org.transept.operator.TmSafeGapShiftMutation;
import org.transept.structure.ProteinStructures;
import org.uma.evolver.parameter.ConditionalParameter;
import org.uma.evolver.parameter.catalogue.mutationparameter.MutationParameter;
import org.uma.jmetal.operator.mutation.MutationOperator;

/**
 * Mutation of a parameter space, with the global sub-parameter {@code mutationProbability} and,
 * for the mutations that move gap runs at random, the conditional {@code maxShift}. The algorithm
 * gives the structures and the source of random numbers as the non-configurable sub-parameters
 * {@code structures} and {@code random}.
 */
public class MsaMutationParameter extends MutationParameter<MsaSolution> {

  /**
   * Creates the parameter.
   *
   * @param mutations allowed mutation names
   * @throws IllegalArgumentException if a name is unknown
   */
  public MsaMutationParameter(List<String> mutations) {
    super(mutations);
    MsaParameters.requireKnown(mutations, MsaOperators.MUTATIONS, "mutation");
  }

  @Override
  public MutationOperator<MsaSolution> getMutation() {
    return mutation((Double) findGlobalSubParameter("mutationProbability").value());
  }

  /** Returns the limit of the chosen mutation, if it moves gap runs at random. */
  private int maxShift() {
    return conditionalParameters().stream()
        .filter(conditional -> conditional.condition().test(value()))
        .map(ConditionalParameter::parameter)
        .filter(parameter -> parameter.name().equals("maxShift"))
        .findFirst()
        .map(parameter -> (Integer) parameter.value())
        .orElse(TmSafeGapShiftMutation.UNLIMITED_SHIFT);
  }

  /**
   * Builds the chosen mutation with another probability, e.g. 1 to perturb copies of seeds.
   *
   * @param probability mutation probability
   * @return the mutation
   */
  public MutationOperator<MsaSolution> mutation(double probability) {
    return MsaOperators.mutation(
        value(),
        probability,
        maxShift(),
        (ProteinStructures) nonConfigurableSubParameters().get("structures"),
        (RandomSource) nonConfigurableSubParameters().get("random"));
  }
}
