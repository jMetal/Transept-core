package org.tmm2structalign.algorithm;

/** Mutation operators available to the algorithms. */
public enum MutationType {
  /**
   * The original gap-shift mutation ({@code TmAwareGapShiftMutation}), which only discourages
   * moving gaps into transmembrane helices. Default, to keep the published behaviour.
   */
  GAP_SHIFT,
  /**
   * Gap-shift mutation that never breaks a transmembrane helix ({@code TmSafeGapShiftMutation}),
   * with the same boundary rule as the crossover.
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
  RANDOM_CHOICE;

  /** Tells whether the operator can break a transmembrane helix. */
  public boolean isTmSafe() {
    return this != GAP_SHIFT;
  }
}
