package org.transept.evolver;

import org.transept.alignment.MsaSolution;
import org.uma.evolver.parameter.ParameterSpace;
import org.uma.jmetal.component.catalogue.common.evaluation.Evaluation;
import org.uma.jmetal.component.catalogue.common.evaluation.impl.MultiThreadedEvaluation;
import org.uma.jmetal.component.catalogue.common.evaluation.impl.SequentialEvaluation;
import org.uma.jmetal.component.catalogue.common.evaluation.impl.SequentialEvaluationWithArchive;
import org.uma.jmetal.component.catalogue.common.termination.Termination;
import org.uma.jmetal.component.catalogue.common.termination.impl.TerminationByComputingTime;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.util.archive.Archive;

/** What the population-based MSA algorithms built with Evolver have in common. */
final class MsaComponents {

  private MsaComponents() {}

  /**
   * Gives the operators and the initial population of a parameter space what they cannot be
   * configured with: the source of random numbers, the structures and the seeds.
   */
  static void setNonConfigurableParameters(ParameterSpace space, MsaRunContext context) {
    space.get("crossover").addNonConfigurableSubParameter("random", context.random());
    var mutation = (MsaMutationParameter) space.get("mutation");
    mutation.addNonConfigurableSubParameter("structures", context.structures());
    mutation.addNonConfigurableSubParameter("random", context.random());
    space
        .get("createInitialSolutions")
        .addNonConfigurableSubParameter("seeds", context.seeds())
        .addNonConfigurableSubParameter("random", context.random())
        .addNonConfigurableSubParameter("copyMutation", mutation.mutation(1.0));
  }

  /**
   * Evaluation with as many threads as the context asks for, feeding the archive if there is one
   * (Evolver requires a {@link SequentialEvaluationWithArchive} then).
   */
  static Evaluation<MsaSolution> evaluation(
      Problem<MsaSolution> problem, Archive<MsaSolution> archive, MsaRunContext context) {
    int threads = context.threads();
    if (archive != null) {
      return threads > 1
          ? new ParallelEvaluationWithArchive(threads, problem, archive)
          : new SequentialEvaluationWithArchive<>(problem, archive);
    }
    return threads > 1
        ? new MultiThreadedEvaluation<>(threads, problem)
        : new SequentialEvaluation<>(problem);
  }

  /** The time limit of the context, or the given termination by evaluations. */
  static Termination termination(MsaRunContext context, Termination byEvaluations) {
    return context.maxTime() == null
        ? byEvaluations
        : new TerminationByComputingTime(context.maxTime().toMillis());
  }
}
