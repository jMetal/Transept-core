package org.tmm2structalign.operator;

/**
 * Opens a new gap run between two residues of every row without breaking a transmembrane helix.
 *
 * <p>TM-safe version of {@code SplitANonGapsGroupMSAMutation} (SANGG) of M2Align and
 * TM-MSAligner. The preferred insertion positions lie between two residue columns, so the
 * inserted gaps form a new run instead of extending an existing one. Every such position is
 * equally likely, including those before the first gap run, which the original never chose (it
 * mistook them for positions inside a run). Rows without a valid position between two residues
 * receive the gaps at any valid position. The alignment grows by {@code gapLength} columns.
 */
public final class TmSafeSplitResidueRunMutation extends TmSafeGapInsertion {

  /**
   * Creates the operator with one gap per row and the global jMetal generator.
   *
   * @param probability probability of mutating a solution
   */
  public TmSafeSplitResidueRunMutation(double probability) {
    this(probability, 1, RandomSource.jmetal());
  }

  /**
   * Creates the operator.
   *
   * @param probability probability of mutating a solution
   * @param gapLength gaps inserted into every row
   * @param random source of random numbers
   */
  public TmSafeSplitResidueRunMutation(double probability, int gapLength, RandomSource random) {
    super(probability, gapLength, random);
  }

  @Override
  boolean isPreferred(boolean[] gaps, int position) {
    return position > 0 && position < gaps.length && !gaps[position - 1] && !gaps[position];
  }
}
