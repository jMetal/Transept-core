package org.tmm2structalign.algorithm;

/** Crossover operators available to the algorithms. Both never cut a transmembrane helix. */
public enum CrossoverType {
  /** Single-point crossover ({@code TmSafeSinglePointCrossover}). Default, as published. */
  SINGLE_POINT,
  /**
   * Two-point crossover ({@code TmSafeMultiPointCrossover}): the child takes a middle segment
   * from the second parent.
   */
  TWO_POINT
}
