package org.transept.evolver;

import java.util.List;
import org.transept.alignment.MsaSolution;
import org.transept.operator.RandomSource;
import org.transept.operator.TmSafeBlockShiftMutation;
import org.transept.operator.TmSafeGapShiftMutation;
import org.transept.structure.ProteinStructures;
import org.uma.evolver.parameter.ConditionalParameter;
import org.uma.evolver.parameter.catalogue.mutationparameter.MutationParameter;
import org.uma.jmetal.operator.mutation.MutationOperator;

/**
 * Mutation of a parameter space, with the global sub-parameter {@code mutationProbability} and, for
 * the mutations that have them, the conditional {@code maxShift} and {@code blockShift}. The
 * algorithm gives the structures and the source of random numbers as the non-configurable
 * sub-parameters {@code structures} and {@code random}.
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

  /** Returns a setting of the chosen mutation, or the given value if it does not have it. */
  private int conditionalValue(String name, int otherwise) {
    return conditionalParameters().stream()
        .filter(conditional -> conditional.condition().test(value()))
        .map(ConditionalParameter::parameter)
        .filter(parameter -> parameter.name().equals(name))
        .findFirst()
        .map(parameter -> (Integer) parameter.value())
        .orElse(otherwise);
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
        conditionalValue("maxShift", TmSafeGapShiftMutation.UNLIMITED_SHIFT),
        conditionalValue("blockShift", TmSafeBlockShiftMutation.DEFAULT_MAX_SHIFT),
        (ProteinStructures) nonConfigurableSubParameters().get("structures"),
        (RandomSource) nonConfigurableSubParameters().get("random"));
  }
}
