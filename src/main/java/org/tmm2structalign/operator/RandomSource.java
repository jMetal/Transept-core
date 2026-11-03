package org.tmm2structalign.operator;

import org.uma.jmetal.util.pseudorandom.JMetalRandom;

/** Source of random numbers for the variation operators. */
public interface RandomSource {

  /**
   * Returns a uniformly distributed integer.
   *
   * @param lower lower bound, inclusive
   * @param upper upper bound, inclusive
   * @return an integer in {@code [lower, upper]}
   */
  int nextInt(int lower, int upper);

  /** Returns a uniformly distributed double in {@code [0, 1)}. */
  double nextDouble();

  /**
   * Returns the global jMetal generator, so that runs are reproducible from the seed given to
   * {@link JMetalRandom#setSeed(long)}.
   */
  static RandomSource jmetal() {
    return new RandomSource() {
      @Override
      public int nextInt(int lower, int upper) {
        return JMetalRandom.getInstance().nextInt(lower, upper);
      }

      @Override
      public double nextDouble() {
        return JMetalRandom.getInstance().nextDouble();
      }
    };
  }
}
