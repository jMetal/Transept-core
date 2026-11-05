package org.transept.objective;

/**
 * Affine gap penalties, subtracted from the score.
 *
 * @param open penalty of the first gap position of a run
 * @param extend penalty of every further gap position of the run
 */
public record GapPenalties(double open, double extend) {

  /** Validates that penalties are finite and non-negative. */
  public GapPenalties {
    if (!Double.isFinite(open) || !Double.isFinite(extend) || open < 0 || extend < 0) {
      throw new IllegalArgumentException(
          "Gap penalties must be finite and non-negative: open=" + open + ", extend=" + extend);
    }
  }
}
