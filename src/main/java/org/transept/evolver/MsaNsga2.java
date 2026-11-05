package org.transept.evolver;

import org.transept.alignment.MsaSolution;
import org.uma.evolver.algorithm.BaseLevelAlgorithm;
import org.uma.evolver.algorithm.nsgaii.BaseNSGAII;
import org.uma.evolver.parameter.ParameterSpace;
import org.uma.jmetal.component.catalogue.common.evaluation.Evaluation;
import org.uma.jmetal.component.catalogue.common.termination.Termination;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.util.archive.Archive;

/**
 * NSGA-II for multiple sequence alignment, configured from the parameter space
 * {@code MSANSGAII.yaml} (see {@link MsaConfigurations}). Compared with Evolver's NSGA-II for
 * permutations, it starts from the seed alignments, evaluates with several threads and can stop
 * after a time limit; the rest (selection, replacement, external archive) is Evolver's.
 */
public class MsaNsga2 extends BaseNSGAII<MsaSolution> {

  /** Name of the parameter space resource. */
  public static final String PARAMETER_SPACE = "MSANSGAII.yaml";

  private final MsaRunContext context;

  /**
   * Creates an algorithm template; {@link #createInstance} gives runnable instances.
   *
   * @param populationSize population size
   * @param parameterSpace parameter space of the algorithm
   * @param context seeds, structures, random numbers, threads and time limit
   */
  public MsaNsga2(int populationSize, ParameterSpace parameterSpace, MsaRunContext context) {
    super(populationSize, parameterSpace);
    this.context = context;
  }

  /**
   * Creates a runnable algorithm.
   *
   * @param problem problem to solve
   * @param populationSize population size
   * @param maximumNumberOfEvaluations stopping criterion, unless the context sets a time limit
   * @param parameterSpace parameter space of the algorithm
   * @param context seeds, structures, random numbers, threads and time limit
   */
  public MsaNsga2(
      Problem<MsaSolution> problem,
      int populationSize,
      int maximumNumberOfEvaluations,
      ParameterSpace parameterSpace,
      MsaRunContext context) {
    super(problem, populationSize, maximumNumberOfEvaluations, parameterSpace);
    this.context = context;
  }

  @Override
  public BaseLevelAlgorithm<MsaSolution> createInstance(
      Problem<MsaSolution> problem, int maximumNumberOfEvaluations) {
    return new MsaNsga2(
        problem,
        populationSize,
        maximumNumberOfEvaluations,
        parameterSpace.createInstance(),
        context);
  }

  @Override
  protected void setNonConfigurableParameters() {
    MsaComponents.setNonConfigurableParameters(parameterSpace, context);
  }

  @Override
  protected Evaluation<MsaSolution> createEvaluation(Archive<MsaSolution> archive) {
    return MsaComponents.evaluation(problem, archive, context);
  }

  @Override
  protected Termination createTermination() {
    return MsaComponents.termination(context, super.createTermination());
  }
}
