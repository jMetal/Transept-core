package org.transept.evolver;

import org.transept.alignment.MsaSolution;
import org.uma.evolver.algorithm.BaseLevelAlgorithm;
import org.uma.evolver.algorithm.rdemoea.BaseRDEMOEA;
import org.uma.evolver.parameter.ParameterSpace;
import org.uma.jmetal.component.catalogue.common.evaluation.Evaluation;
import org.uma.jmetal.component.catalogue.common.termination.Termination;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.util.archive.Archive;

/**
 * SPEA2-style algorithm for multiple sequence alignment, configured from the parameter space
 * {@code MSASPEA2.yaml}. It is Evolver's RDEMOEA (a generational algorithm with a configurable
 * ranking and density estimator) whose default configuration is that of SPEA2: strength ranking,
 * k-nearest-neighbour density and sequential truncation of the worst solutions. Unlike the
 * original SPEA2, the population plays the role of the archive. Choosing dominance ranking and
 * crowding distance turns it into NSGA-II.
 */
public class MsaSpea2 extends BaseRDEMOEA<MsaSolution> {

  /** Name of the parameter space resource. */
  public static final String PARAMETER_SPACE = "MSASPEA2.yaml";

  private final MsaRunContext context;

  /**
   * Creates an algorithm template; {@link #createInstance} gives runnable instances.
   *
   * @param populationSize population size
   * @param parameterSpace parameter space of the algorithm
   * @param context seeds, structures, random numbers, threads and time limit
   */
  public MsaSpea2(int populationSize, ParameterSpace parameterSpace, MsaRunContext context) {
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
  public MsaSpea2(
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
    return new MsaSpea2(
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
