package org.transept.evolver;

import java.util.Comparator;
import org.transept.alignment.MsaSolution;
import org.uma.evolver.algorithm.BaseLevelAlgorithm;
import org.uma.evolver.parameter.ParameterSpace;
import org.uma.evolver.parameter.catalogue.createinitialsolutionsparameter.CreateInitialSolutionsParameter;
import org.uma.evolver.parameter.catalogue.selectionparameter.SelectionParameter;
import org.uma.evolver.parameter.catalogue.variationparameter.VariationParameter;
import org.uma.jmetal.component.algorithm.EvolutionaryAlgorithm;
import org.uma.jmetal.component.catalogue.common.termination.impl.TerminationByEvaluations;
import org.uma.jmetal.component.catalogue.ea.replacement.Replacement;
import org.uma.jmetal.component.catalogue.ea.replacement.impl.MuCommaLambdaReplacement;
import org.uma.jmetal.component.catalogue.ea.replacement.impl.MuPlusLambdaReplacement;
import org.uma.jmetal.component.catalogue.ea.variation.Variation;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.util.comparator.ObjectiveComparator;

/**
 * Generational genetic algorithm for the single-objective formulation (the lDDT alone),
 * configured from the parameter space {@code MSAGA.yaml}: the seeded initial population and the
 * operators of the multiobjective algorithms, a selection on the objective and a (mu + lambda) or
 * (mu, lambda) replacement ({@code gaReplacement}). Evolver has no base-level genetic algorithm,
 * so this class assembles it from Evolver's parameters and jMetal's components.
 */
public class MsaGeneticAlgorithm implements BaseLevelAlgorithm<MsaSolution> {

  /** Name of the parameter space resource. */
  public static final String PARAMETER_SPACE = "MSAGA.yaml";

  /** Name reported by jMetal. */
  public static final String NAME = "GeneticAlgorithm";

  private final ParameterSpace parameterSpace;
  private final MsaRunContext context;
  private final int populationSize;
  private Problem<MsaSolution> problem;
  private int maximumNumberOfEvaluations;

  /**
   * Creates an algorithm template; {@link #createInstance} gives runnable instances.
   *
   * @param populationSize population size
   * @param parameterSpace parameter space of the algorithm
   * @param context seeds, structures, random numbers, threads and time limit
   */
  public MsaGeneticAlgorithm(
      int populationSize, ParameterSpace parameterSpace, MsaRunContext context) {
    this.populationSize = populationSize;
    this.parameterSpace = parameterSpace;
    this.context = context;
  }

  /**
   * Creates a runnable algorithm.
   *
   * @param problem problem to solve; it must have a single objective
   * @param populationSize population size
   * @param maximumNumberOfEvaluations stopping criterion, unless the context sets a time limit
   * @param parameterSpace parameter space of the algorithm
   * @param context seeds, structures, random numbers, threads and time limit
   */
  public MsaGeneticAlgorithm(
      Problem<MsaSolution> problem,
      int populationSize,
      int maximumNumberOfEvaluations,
      ParameterSpace parameterSpace,
      MsaRunContext context) {
    this(populationSize, parameterSpace, context);
    this.problem = problem;
    this.maximumNumberOfEvaluations = maximumNumberOfEvaluations;
  }

  @Override
  public ParameterSpace parameterSpace() {
    return parameterSpace;
  }

  @Override
  public BaseLevelAlgorithm<MsaSolution> createInstance(
      Problem<MsaSolution> problem, int maximumNumberOfEvaluations) {
    return new MsaGeneticAlgorithm(
        problem,
        populationSize,
        maximumNumberOfEvaluations,
        parameterSpace.createInstance(),
        context);
  }

  @Override
  @SuppressWarnings("unchecked")
  public EvolutionaryAlgorithm<MsaSolution> build() {
    if (problem.numberOfObjectives() != 1) {
      throw new IllegalArgumentException(
          "The genetic algorithm optimises a single objective, got "
              + problem.numberOfObjectives());
    }
    MsaComponents.setNonConfigurableParameters(parameterSpace, context);
    var variationParameter = (VariationParameter<MsaSolution>) parameterSpace.get("variation");
    variationParameter.addNonConfigurableSubParameter(
        "offspringPopulationSize", parameterSpace.get("offspringPopulationSize").value());
    Variation<MsaSolution> variation = variationParameter.getVariation();
    Comparator<MsaSolution> fitness = new ObjectiveComparator<>(0);
    var selection =
        ((SelectionParameter<MsaSolution>) parameterSpace.get("selection"))
            .getSelection(variation.matingPoolSize(), fitness);
    var creationParameter =
        (CreateInitialSolutionsParameter<MsaSolution>) parameterSpace.get("createInitialSolutions");
    var creation = creationParameter.getCreateInitialSolutionsStrategy(problem, populationSize);
    return new EvolutionaryAlgorithm<>(
        NAME,
        creation,
        MsaComponents.evaluation(problem, null, context),
        MsaComponents.termination(
            context, new TerminationByEvaluations(maximumNumberOfEvaluations)),
        selection,
        variation,
        replacement(fitness));
  }

  private Replacement<MsaSolution> replacement(Comparator<MsaSolution> fitness) {
    String value = (String) parameterSpace.get("gaReplacement").value();
    return switch (value) {
      case "muPlusLambda" -> new MuPlusLambdaReplacement<>(fitness);
      case "muCommaLambda" -> new MuCommaLambdaReplacement<>(fitness);
      default -> throw new IllegalArgumentException("Unknown gaReplacement: " + value);
    };
  }
}
