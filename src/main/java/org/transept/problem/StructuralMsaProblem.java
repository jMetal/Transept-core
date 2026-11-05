package org.transept.problem;

import java.util.List;
import org.transept.alignment.Msa;
import org.transept.alignment.MsaSolution;
import org.transept.model.ProteinSet;
import org.transept.objective.Objective;
import org.uma.jmetal.problem.Problem;

/**
 * Multiple sequence alignment of a fixed protein set, evaluated with a list of objectives.
 *
 * <p>Evaluation first removes gap-only columns from the solution (it modifies the solution) and
 * then stores every objective value in jMetal's minimisation convention: maximised objectives are
 * negated.
 *
 * <p>Solutions are created by a seeded population (see
 * {@code org.transept.algorithm.SeededPopulationCreation}), so
 * {@link #createSolution()} is not supported.
 */
public final class StructuralMsaProblem implements Problem<MsaSolution> {

  private static final long serialVersionUID = 1L;

  private final String name;
  private final transient ProteinSet proteins;
  private final transient List<Objective> objectives;

  /**
   * Creates the problem.
   *
   * @param name problem name, usually the dataset name
   * @param proteins proteins to align
   * @param objectives objectives in output order; at least one
   */
  public StructuralMsaProblem(String name, ProteinSet proteins, List<Objective> objectives) {
    if (objectives == null || objectives.isEmpty()) {
      throw new IllegalArgumentException("At least one objective is required");
    }
    this.name = name;
    this.proteins = proteins;
    this.objectives = List.copyOf(objectives);
  }

  /** Returns the proteins to align. */
  public ProteinSet proteins() {
    return proteins;
  }

  /** Returns the objectives in output order. */
  public List<Objective> objectives() {
    return objectives;
  }

  @Override
  public int numberOfVariables() {
    return proteins.size();
  }

  @Override
  public int numberOfObjectives() {
    return objectives.size();
  }

  @Override
  public int numberOfConstraints() {
    return 0;
  }

  @Override
  public String name() {
    return name;
  }

  @Override
  public MsaSolution evaluate(MsaSolution solution) {
    solution.compactGapColumns();
    Msa msa = solution.toMsa();
    for (int index = 0; index < objectives.size(); index++) {
      Objective objective = objectives.get(index);
      double value = objective.evaluate(msa);
      solution.objectives()[index] = objective.isMaximized() ? -value : value;
    }
    return solution;
  }

  /**
   * Converts a stored objective value back to its natural direction.
   *
   * @param solution evaluated solution
   * @param index objective index
   * @return the objective value as reported to users
   */
  public double naturalValue(MsaSolution solution, int index) {
    double stored = solution.objectives()[index];
    return objectives.get(index).isMaximized() ? -stored : stored;
  }

  /**
   * Creates a solution from an alignment of this problem's proteins.
   *
   * @param msa alignment whose protein set must be this problem's
   * @return a new, unevaluated solution
   */
  public MsaSolution createSolution(Msa msa) {
    if (msa.proteins() != proteins) {
      throw new IllegalArgumentException("The alignment belongs to a different protein set");
    }
    return new MsaSolution(proteins, msa.rows(), objectives.size());
  }

  /**
   * Not supported: solutions are created from seed alignments.
   *
   * @throws UnsupportedOperationException always
   */
  @Override
  public MsaSolution createSolution() {
    throw new UnsupportedOperationException(
        "Random solutions are not supported; create solutions from seed alignments");
  }
}
