package org.tmm2structalign.operator;

/**
 * Inserts gaps at a random position of every row without breaking a transmembrane helix.
 *
 * <p>TM-safe version of {@code InsertARandomGapMSAMutation} (IRG) of M2Align and TM-MSAligner:
 * every insertion position of a row that does not separate two transmembrane residues is equally
 * likely, including those that extend an existing gap run. The alignment grows by
 * {@code gapLength} columns.
 */
public final class TmSafeInsertGapMutation extends TmSafeGapInsertion {

  /**
   * Creates the operator with one gap per row and the global jMetal generator.
   *
   * @param probability probability of mutating a solution
   */
  public TmSafeInsertGapMutation(double probability) {
    this(probability, 1, RandomSource.jmetal());
  }

  /**
   * Creates the operator.
   *
   * @param probability probability of mutating a solution
   * @param gapLength gaps inserted into every row
   * @param random source of random numbers
   */
  public TmSafeInsertGapMutation(double probability, int gapLength, RandomSource random) {
    super(probability, gapLength, random);
  }

  @Override
  boolean isPreferred(boolean[] gaps, int position) {
    return true;
  }
}
