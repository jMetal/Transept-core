package org.transept.evolver;

import java.util.List;
import org.transept.algorithm.SeededPopulationCreation;
import org.transept.alignment.Msa;
import org.transept.alignment.MsaSolution;
import org.transept.operator.RandomSource;
import org.transept.problem.StructuralMsaProblem;
import org.uma.evolver.parameter.catalogue.createinitialsolutionsparameter.CreateInitialSolutionsParameter;
import org.uma.jmetal.component.catalogue.common.solutionscreation.SolutionsCreation;
import org.uma.jmetal.operator.mutation.MutationOperator;
import org.uma.jmetal.problem.Problem;

/**
 * Initial population of a parameter space: {@code seededDistinct}, the seed alignments followed
 * by distinct crossover children ({@link SeededPopulationCreation}). The algorithm gives the
 * seeds, the source of random numbers and the mutation for the fallback copies as the
 * non-configurable sub-parameters {@code seeds}, {@code random} and {@code copyMutation}.
 */
public class CreateInitialSolutionsMsaParameter
    extends CreateInitialSolutionsParameter<MsaSolution> {

  /** Seeds followed by distinct crossover children. */
  public static final String SEEDED_DISTINCT = "seededDistinct";

  /**
   * Creates the parameter.
   *
   * @param strategies allowed strategies
   * @throws IllegalArgumentException if a strategy is unknown
   */
  public CreateInitialSolutionsMsaParameter(List<String> strategies) {
    super(DEFAULT_NAME, strategies);
    MsaParameters.requireKnown(strategies, List.of(SEEDED_DISTINCT), DEFAULT_NAME);
  }

  @Override
  @SuppressWarnings("unchecked")
  public SolutionsCreation<MsaSolution> getCreateInitialSolutionsStrategy(
      Problem<MsaSolution> problem, int populationSize) {
    return new SeededPopulationCreation(
        (StructuralMsaProblem) problem,
        (List<Msa>) nonConfigurableSubParameters().get("seeds"),
        populationSize,
        (RandomSource) nonConfigurableSubParameters().get("random"),
        (MutationOperator<MsaSolution>) nonConfigurableSubParameters().get("copyMutation"),
        (Integer) nonConfigurableSubParameters().get("cutCandidates"));
  }
}
