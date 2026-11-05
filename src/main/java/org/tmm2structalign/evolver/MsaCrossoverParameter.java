package org.tmm2structalign.evolver;

import java.util.List;
import org.tmm2structalign.alignment.MsaSolution;
import org.tmm2structalign.operator.RandomSource;
import org.uma.evolver.parameter.catalogue.crossoverparameter.CrossoverParameter;
import org.uma.jmetal.operator.crossover.CrossoverOperator;

/**
 * Crossover of a parameter space, with the global sub-parameter {@code crossoverProbability}.
 * The algorithm gives the source of random numbers as the non-configurable sub-parameter
 * {@code random}.
 */
public class MsaCrossoverParameter extends CrossoverParameter<MsaSolution> {

  /**
   * Creates the parameter.
   *
   * @param crossovers allowed crossover names
   * @throws IllegalArgumentException if a name is unknown
   */
  public MsaCrossoverParameter(List<String> crossovers) {
    super(crossovers);
    MsaParameters.requireKnown(crossovers, MsaOperators.CROSSOVERS, "crossover");
  }

  @Override
  public CrossoverOperator<MsaSolution> getCrossover() {
    double probability = (Double) findGlobalSubParameter("crossoverProbability").value();
    return MsaOperators.crossover(
        value(), probability, (RandomSource) nonConfigurableSubParameters().get("random"));
  }
}
