package org.transept.algorithm;

/** Crossover operators available to the algorithms. None of them cuts a transmembrane helix. */
public enum CrossoverType {
  /**
   * Single-point crossover that never cuts a transmembrane helix ({@code
   * TmSafeSinglePointCrossover}). Default.
   */
  TM_SAFE_SINGLE_POINT,
  /**
   * Uniform crossover between the positions where both parents agree, which adds no gap ({@code
   * TmSafeUniformCrossover}).
   */
  TM_SAFE_UNIFORM
}
