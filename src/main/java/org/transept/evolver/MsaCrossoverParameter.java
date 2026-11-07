package org.transept.evolver;

import java.util.List;
import org.transept.alignment.MsaSolution;
import org.transept.operator.RandomSource;
import org.uma.evolver.parameter.catalogue.crossoverparameter.CrossoverParameter;
import org.uma.jmetal.operator.crossover.CrossoverOperator;

/**
 * Crossover of a parameter space, with the global sub-parameters {@code crossoverProbability} and
 * {@code cutCandidates}. The algorithm gives the source of random numbers as the non-configurable
 * sub-parameter {@code random}.
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

  /**
   * Returns how many valid cuts the crossover draws to use the one with the least padding.
   *
   * @return the value of the sub-parameter {@code cutCandidates}
   */
  public int cutCandidates() {
    return (Integer) findGlobalSubParameter("cutCandidates").value();
  }

  @Override
  public CrossoverOperator<MsaSolution> getCrossover() {
    double probability = (Double) findGlobalSubParameter("crossoverProbability").value();
    return MsaOperators.crossover(
        value(),
        probability,
        cutCandidates(),
        (RandomSource) nonConfigurableSubParameters().get("random"));
  }
}
