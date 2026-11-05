package org.transept.objective;

import org.transept.alignment.Msa;

/** Quality measure of a multiple sequence alignment. */
public interface Objective {

  /** Returns a short, stable name used in output headers (e.g. {@code SOPwTP}). */
  String name();

  /** Returns true if higher values are better. */
  boolean isMaximized();

  /**
   * Computes the raw value of the objective.
   *
   * @param msa alignment to evaluate; implementations must not assume gap-only columns were
   *     removed
   * @return the objective value, in the natural direction given by {@link #isMaximized()}
   */
  double evaluate(Msa msa);
}
