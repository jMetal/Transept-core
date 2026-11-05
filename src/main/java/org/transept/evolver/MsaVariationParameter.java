package org.transept.evolver;

import java.util.List;
import org.transept.alignment.MsaSolution;
import org.uma.evolver.parameter.catalogue.crossoverparameter.CrossoverParameter;
import org.uma.evolver.parameter.catalogue.mutationparameter.MutationParameter;
import org.uma.evolver.parameter.catalogue.variationparameter.VariationParameter;
import org.uma.jmetal.component.catalogue.ea.variation.Variation;
import org.uma.jmetal.component.catalogue.ea.variation.impl.CrossoverAndMutationVariation;

/**
 * Variation of a parameter space: {@code crossoverAndMutationVariation}, with the conditional
 * parameters {@code crossover} and {@code mutation}. The algorithm gives the number of offspring
 * as the non-configurable sub-parameter {@code offspringPopulationSize}.
 */
public class MsaVariationParameter extends VariationParameter<MsaSolution> {

  /** The only variation strategy. */
  public static final String CROSSOVER_AND_MUTATION = "crossoverAndMutationVariation";

  /**
   * Creates the parameter.
   *
   * @param strategies allowed variation strategies
   * @throws IllegalArgumentException if a strategy is unknown
   */
  public MsaVariationParameter(List<String> strategies) {
    super(strategies);
    MsaParameters.requireKnown(strategies, List.of(CROSSOVER_AND_MUTATION), "variation");
  }

  @Override
  @SuppressWarnings("unchecked")
  public Variation<MsaSolution> getVariation() {
    var crossover = (CrossoverParameter<MsaSolution>) findConditionalParameter("crossover");
    var mutation = (MutationParameter<MsaSolution>) findConditionalParameter("mutation");
    return new CrossoverAndMutationVariation<>(
        (Integer) nonConfigurableSubParameters().get("offspringPopulationSize"),
        crossover.getCrossover(),
        mutation.getMutation());
  }
}
