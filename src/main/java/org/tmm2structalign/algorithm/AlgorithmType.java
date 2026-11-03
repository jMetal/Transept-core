package org.tmm2structalign.algorithm;

/** Search algorithms available to align a dataset. */
public enum AlgorithmType {
  /** NSGA-II with the published configuration (default). */
  NSGA_II(false, true),
  /** NSGA-III; the population size must be a multiple of 4. */
  NSGA_III(false, true),
  /** AGE-MOEA. */
  AGE_MOEA(false, true),
  /** AGE-MOEA-II. */
  AGE_MOEA_II(false, true),
  /** SMS-EMOA, steady state: one offspring per iteration. */
  SMS_EMOA(true, true),
  /** MOEA/D with normalised PBI aggregation, steady state. */
  MOEA_D(true, true),
  /** PAES, a (1+1) evolution strategy with a crowding-distance archive; mutation only. */
  PAES(true, false),
  /** Random search over alignments with gaps placed at random (weak baseline). */
  RANDOM_SEARCH(false, false),
  /** Random search over seeds perturbed with random mutations (strong baseline). */
  RANDOM_SEEDED(false, false);

  private final boolean steadyState;
  private final boolean usesCrossover;

  AlgorithmType(boolean steadyState, boolean usesCrossover) {
    this.steadyState = steadyState;
    this.usesCrossover = usesCrossover;
  }

  /** Returns true if one solution is evaluated per iteration, so threads do not help. */
  public boolean isSteadyState() {
    return steadyState;
  }

  /** Returns true if the algorithm applies the crossover operator. */
  public boolean usesCrossover() {
    return usesCrossover;
  }
}
