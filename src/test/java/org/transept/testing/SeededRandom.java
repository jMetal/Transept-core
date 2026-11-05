package org.transept.testing;

import java.util.Random;
import org.transept.operator.RandomSource;

/** Reproducible random source for tests, independent of the global jMetal generator. */
public final class SeededRandom implements RandomSource {

  private final Random random;

  /**
   * Creates the source.
   *
   * @param seed seed of the underlying generator
   */
  public SeededRandom(long seed) {
    this.random = new Random(seed);
  }

  @Override
  public int nextInt(int lower, int upper) {
    return lower + random.nextInt(upper - lower + 1);
  }

  @Override
  public double nextDouble() {
    return random.nextDouble();
  }
}
