package org.transept.alignment;

import java.util.List;
import org.transept.model.ProteinSet;
import org.uma.jmetal.solution.AbstractSolution;

/**
 * A jMetal solution holding one gap layout per protein, in canonical order.
 *
 * <p>Variables are immutable values: operators replace them with {@code variables().set(...)}
 * instead of mutating them.
 */
public class MsaSolution extends AbstractSolution<GapIntervals> {

  private static final long serialVersionUID = 1L;

  private final transient ProteinSet proteins;

  /**
   * Creates a solution from one gap layout per protein.
   *
   * @param proteins proteins in canonical order
   * @param rows gap layouts, validated against the proteins
   * @param numberOfObjectives number of objectives of the problem
   */
  public MsaSolution(ProteinSet proteins, List<GapIntervals> rows, int numberOfObjectives) {
    super(proteins.size(), numberOfObjectives);
    Msa.of(proteins, rows);
    this.proteins = proteins;
    for (int sequence = 0; sequence < rows.size(); sequence++) {
      variables().set(sequence, rows.get(sequence));
    }
  }

  private MsaSolution(MsaSolution solution) {
    super(solution.variables().size(), solution.objectives().length);
    this.proteins = solution.proteins;
    for (int sequence = 0; sequence < solution.variables().size(); sequence++) {
      variables().set(sequence, solution.variables().get(sequence));
    }
    System.arraycopy(
        solution.objectives(), 0, objectives(), 0, solution.objectives().length);
    attributes.putAll(solution.attributes());
  }

  @Override
  public MsaSolution copy() {
    return new MsaSolution(this);
  }

  /** Returns the proteins in canonical order. */
  public ProteinSet proteins() {
    return proteins;
  }

  /** Returns the number of columns of the alignment. */
  public int alignmentLength() {
    return variables().get(0).alignedLength();
  }

  /**
   * Decodes the solution into an alignment view.
   *
   * @return the alignment
   * @throws IllegalArgumentException if the variables were replaced by inconsistent rows
   */
  public Msa toMsa() {
    return Msa.of(proteins, variables());
  }

  /** Removes the columns in which every row has a gap. */
  public void compactGapColumns() {
    List<GapIntervals> compacted = GapColumnCompactor.compact(variables());
    for (int sequence = 0; sequence < compacted.size(); sequence++) {
      variables().set(sequence, compacted.get(sequence));
    }
  }
}
