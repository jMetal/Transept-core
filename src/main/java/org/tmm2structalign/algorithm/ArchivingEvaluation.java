package org.tmm2structalign.algorithm;

import java.util.List;
import org.tmm2structalign.alignment.MsaSolution;
import org.uma.jmetal.component.catalogue.common.evaluation.Evaluation;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.util.archive.Archive;

/**
 * Evaluation that also offers a copy of every evaluated solution to an external archive, so that
 * the result of a run can be the archive instead of the final population.
 */
public final class ArchivingEvaluation implements Evaluation<MsaSolution> {

  private final Evaluation<MsaSolution> evaluation;
  private final Archive<MsaSolution> archive;

  /**
   * Creates the evaluation.
   *
   * @param evaluation evaluation that computes the objectives
   * @param archive archive that receives a copy of every evaluated solution
   */
  public ArchivingEvaluation(Evaluation<MsaSolution> evaluation, Archive<MsaSolution> archive) {
    this.evaluation = evaluation;
    this.archive = archive;
  }

  @Override
  public List<MsaSolution> evaluate(List<MsaSolution> solutions) {
    List<MsaSolution> evaluated = evaluation.evaluate(solutions);
    evaluated.forEach(solution -> archive.add(solution.copy()));
    return evaluated;
  }

  @Override
  public int computedEvaluations() {
    return evaluation.computedEvaluations();
  }

  @Override
  public Problem<MsaSolution> problem() {
    return evaluation.problem();
  }

  /** Returns the archive. */
  public Archive<MsaSolution> archive() {
    return archive;
  }
}
