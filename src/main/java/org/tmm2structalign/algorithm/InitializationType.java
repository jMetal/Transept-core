package org.tmm2structalign.algorithm;

/** How the non-seed members of the initial population are created. */
public enum InitializationType {
  /**
   * Children of the TM-safe crossover between random members, as in the original
   * implementation. Default, to keep the published behaviour. Duplicates are possible.
   */
  CROSSOVER,
  /**
   * Like {@link #CROSSOVER}, but a child is kept only if it differs from every member once
   * gap-only columns are removed. If not enough distinct children are found, the remaining
   * members are perturbed with mutations.
   */
  CROSSOVER_DISTINCT,
  /** Children of the crossover followed by several mutations, for more local diversity. */
  PERTURBED
}
