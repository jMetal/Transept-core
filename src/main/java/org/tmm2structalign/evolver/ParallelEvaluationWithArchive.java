package org.tmm2structalign.evolver;

import java.util.List;
import java.util.concurrent.ForkJoinPool;
import org.tmm2structalign.alignment.MsaSolution;
import org.uma.jmetal.component.catalogue.common.evaluation.impl.SequentialEvaluationWithArchive;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.util.archive.Archive;

/**
 * Evaluates the solutions with several threads and offers a copy of each to an archive, in the
 * order of the list, so the archive does not depend on the number of threads. It extends
 * {@link SequentialEvaluationWithArchive} because Evolver requires that class when an algorithm
 * returns an external archive.
 */
public class ParallelEvaluationWithArchive extends SequentialEvaluationWithArchive<MsaSolution> {

  private final ForkJoinPool pool;
  private int computedEvaluations;

  /**
   * Creates the evaluation.
   *
   * @param threads number of evaluation threads
   * @param problem problem whose objectives are computed
   * @param archive archive that receives a copy of every evaluated solution
   */
  public ParallelEvaluationWithArchive(
      int threads, Problem<MsaSolution> problem, Archive<MsaSolution> archive) {
    super(problem, archive);
    this.pool = new ForkJoinPool(threads);
  }

  @Override
  public List<MsaSolution> evaluate(List<MsaSolution> solutions) {
    pool.submit(() -> solutions.parallelStream().forEach(problem()::evaluate)).join();
    computedEvaluations = solutions.size();
    solutions.forEach(solution -> archive().add(solution.copy()));
    return solutions;
  }

  @Override
  public int computedEvaluations() {
    return computedEvaluations;
  }
}
