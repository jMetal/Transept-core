package org.transept.evolver;

import java.util.List;
import org.transept.alignment.MsaSolution;
import org.transept.operator.RandomSource;
import org.transept.structure.ProteinStructures;
import org.uma.evolver.parameter.catalogue.mutationparameter.MutationParameter;
import org.uma.jmetal.operator.mutation.MutationOperator;

/**
 * Mutation of a parameter space, with the global sub-parameter {@code mutationProbability}. The
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
        (ProteinStructures) nonConfigurableSubParameters().get("structures"),
        (RandomSource) nonConfigurableSubParameters().get("random"));
  }
}
