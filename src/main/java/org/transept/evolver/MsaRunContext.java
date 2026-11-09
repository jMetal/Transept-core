package org.transept.evolver;

import java.time.Duration;
import java.util.List;
import org.transept.algorithm.StopRequest;
import org.transept.alignment.Msa;
import org.transept.operator.RandomSource;
import org.transept.structure.ProteinStructures;

/**
 * What an MSA algorithm needs besides its configuration: the seed alignments of the dataset, the
 * structures, the source of random numbers, the evaluation threads and, optionally, a time limit.
 *
 * @param seeds seed alignments, in order of preference
 * @param structures structures of the proteins, used by the mutations
 * @param random source of random numbers of the operators and the initial population
 * @param threads evaluation threads; 1 evaluates sequentially
 * @param maxTime time limit of the search, or {@code null} to stop by evaluations
 * @param stop request to stop the search early, keeping its results
 */
public record MsaRunContext(
    List<Msa> seeds,
    ProteinStructures structures,
    RandomSource random,
    int threads,
    Duration maxTime,
    StopRequest stop) {

  /** Validates the context and freezes the seeds. */
  public MsaRunContext {
    if (seeds == null || seeds.isEmpty() || structures == null || random == null) {
      throw new IllegalArgumentException("seeds, structures and random are required");
    }
    if (threads < 1) {
      throw new IllegalArgumentException("threads must be positive: " + threads);
    }
    if (maxTime != null && (maxTime.isZero() || maxTime.isNegative())) {
      throw new IllegalArgumentException("maxTime must be positive: " + maxTime);
    }
    seeds = List.copyOf(seeds);
    stop = stop == null ? StopRequest.NEVER : stop;
  }

  /**
   * Creates a context that cannot be stopped early.
   *
   * @param seeds seed alignments, in order of preference
   * @param structures structures of the proteins, used by the mutations
   * @param random source of random numbers of the operators and the initial population
   * @param threads evaluation threads; 1 evaluates sequentially
   * @param maxTime time limit of the search, or {@code null} to stop by evaluations
   */
  public MsaRunContext(
      List<Msa> seeds,
      ProteinStructures structures,
      RandomSource random,
      int threads,
      Duration maxTime) {
    this(seeds, structures, random, threads, maxTime, StopRequest.NEVER);
  }
}
