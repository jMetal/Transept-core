package org.tmm2structalign.objective;

/** How per-anchor lDDT scores are combined into one alignment score. */
public enum LddtAggregation {
  /**
   * Original definition: for every column, the harmonic mean of its mean structural score and
   * its occupancy (residues / sequences), summed over columns and divided by the alignment
   * length. Penalises long alignments with sparse columns (see docs/review, finding A3).
   */
  HARMONIC_OVER_LENGTH,
  /** Mean of all directional anchor scores; independent of column occupancy and length. */
  PAIR_MEAN
}
