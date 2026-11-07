package org.transept.algorithm;

/** Mutation operators available to the algorithms. None of them breaks a transmembrane helix. */
public enum MutationType {
  /**
   * Gap-shift mutation that never breaks a transmembrane helix ({@code TmSafeGapShiftMutation}),
   * with the same boundary rule as the crossover. Default.
   */
  TM_SAFE_GAP_SHIFT,
  /**
   * Merges two consecutive gap runs of one row ({@code MergeAdjacentGapRunsMutation}, MAGG of
   * TM-MSAligner).
   */
  MERGE_GAP_RUNS,
  /**
   * Applies {@code TM_SAFE_GAP_SHIFT} or {@code MERGE_GAP_RUNS}, chosen uniformly at random each
   * time ({@code RandomChoiceMutation}).
   */
  RANDOM_CHOICE,
  /**
   * Moves a block of gaps a few columns in several rows at once ({@code
   * TmSafeBlockShiftMutation}).
   */
  TM_SAFE_BLOCK_SHIFT
}
