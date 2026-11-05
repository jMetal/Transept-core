package org.transept.objective;

import org.transept.alignment.Msa;
import org.transept.model.ProteinSet;

/**
 * Fraction of residues that the alignment pairs with residues of other proteins; maximised.
 *
 * <p>For every pair of proteins, the number of columns where both have a residue is divided by
 * the length of the shorter protein, the largest number of pairs the two can share. The objective
 * is the sum of the aligned pairs over all pairs of proteins divided by the sum of those maxima,
 * so it lies in [0, 1] and does not need a reference alignment. It counts how much the alignment
 * aligns, not whether it aligns the right residues: the mean of the lDDT anchor scores ignores
 * this, and an alignment that spreads the proteins apart with gaps scores low here.
 */
public final class AlignedPairCoverage implements Objective {

  /** Name used in output files. */
  public static final String NAME = "coverage";

  private final ProteinSet proteins;
  private final double maximumPairs;

  /**
   * Creates the objective for a set of proteins.
   *
   * @param proteins proteins of the alignments that will be evaluated
   * @throws IllegalArgumentException if there are fewer than two proteins
   */
  public AlignedPairCoverage(ProteinSet proteins) {
    if (proteins.size() < 2) {
      throw new IllegalArgumentException("Coverage needs at least two proteins");
    }
    this.proteins = proteins;
    long maximum = 0;
    for (int first = 0; first < proteins.size(); first++) {
      for (int second = first + 1; second < proteins.size(); second++) {
        maximum += Math.min(proteins.get(first).length(), proteins.get(second).length());
      }
    }
    this.maximumPairs = maximum;
  }

  @Override
  public String name() {
    return NAME;
  }

  @Override
  public boolean isMaximized() {
    return true;
  }

  @Override
  public double evaluate(Msa msa) {
    if (msa.proteins() != proteins) {
      throw new IllegalArgumentException("The alignment belongs to a different protein set");
    }
    long aligned = 0;
    for (int column = 0; column < msa.length(); column++) {
      long residues = msa.occupancy(column);
      aligned += residues * (residues - 1) / 2;
    }
    return maximumPairs == 0 ? 0.0 : aligned / maximumPairs;
  }
}
