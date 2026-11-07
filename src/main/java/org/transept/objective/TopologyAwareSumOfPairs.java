package org.transept.objective;

import org.transept.alignment.Msa;

/**
 * Topology-aware sum of pairs (SOPwTP), maximised.
 *
 * <p>For every pair of sequences, residue-residue columns score with the transmembrane matrix
 * when both residues are transmembrane and with the general matrix otherwise. Gaps receive affine
 * penalties, using the transmembrane penalties when either side of the pair is transmembrane
 * (see {@link Msa#isTransmembrane(int, int)} for gaps). Gap-gap columns are skipped and do not
 * close an open gap. Gaps in each sequence of the pair keep separate affine states.
 *
 * <p>Port of {@code StructSumOfPairsWithTopologyPredict}; it reproduces its values exactly.
 */
public final class TopologyAwareSumOfPairs implements Objective {

  /** Name used in output files. */
  public static final String NAME = "SOPwTP";

  private static final int NO_GAP = 0;
  private static final int GAP_IN_FIRST = 1;
  private static final int GAP_IN_SECOND = 2;

  /** Gap penalties of the published experiments when the gap or its partner is transmembrane. */
  public static final GapPenalties TRANSMEMBRANE_GAPS = new GapPenalties(8.0, 3.0);

  /** Gap penalties of the published experiments otherwise. */
  public static final GapPenalties GENERAL_GAPS = new GapPenalties(3.0, 1.0);

  private final SubstitutionMatrix transmembraneMatrix;
  private final SubstitutionMatrix generalMatrix;
  private final GapPenalties transmembraneGaps;
  private final GapPenalties generalGaps;

  /**
   * Creates the objective.
   *
   * @param transmembraneMatrix matrix for transmembrane-transmembrane residue pairs
   * @param generalMatrix matrix for every other residue pair
   * @param transmembraneGaps penalties when the gap or its partner is transmembrane
   * @param generalGaps penalties otherwise
   */
  public TopologyAwareSumOfPairs(
      SubstitutionMatrix transmembraneMatrix,
      SubstitutionMatrix generalMatrix,
      GapPenalties transmembraneGaps,
      GapPenalties generalGaps) {
    this.transmembraneMatrix = transmembraneMatrix;
    this.generalMatrix = generalMatrix;
    this.transmembraneGaps = transmembraneGaps;
    this.generalGaps = generalGaps;
  }

  /** Returns the configuration used in the published experiments (PHAT, BLOSUM62, 8/3, 3/1). */
  public static TopologyAwareSumOfPairs withDefaultParameters() {
    return new TopologyAwareSumOfPairs(
        SubstitutionMatrix.PHAT,
        SubstitutionMatrix.BLOSUM62,
        TRANSMEMBRANE_GAPS,
        GENERAL_GAPS);
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
    double score = 0.0;
    for (int first = 0; first < msa.size() - 1; first++) {
      for (int second = first + 1; second < msa.size(); second++) {
        score += scorePair(msa, first, second);
      }
    }
    return score;
  }

  private double scorePair(Msa msa, int first, int second) {
    double score = 0.0;
    int previousGapState = NO_GAP;
    for (int column = 0; column < msa.length(); column++) {
      boolean gapInFirst = msa.isGap(first, column);
      boolean gapInSecond = msa.isGap(second, column);
      if (gapInFirst && gapInSecond) {
        continue;
      }
      boolean firstTm = msa.isTransmembrane(first, column);
      boolean secondTm = msa.isTransmembrane(second, column);
      if (!gapInFirst && !gapInSecond) {
        SubstitutionMatrix matrix = firstTm && secondTm ? transmembraneMatrix : generalMatrix;
        score += matrix.score(msa.letter(first, column), msa.letter(second, column));
        previousGapState = NO_GAP;
      } else {
        int gapState = gapInFirst ? GAP_IN_FIRST : GAP_IN_SECOND;
        GapPenalties penalties = firstTm || secondTm ? transmembraneGaps : generalGaps;
        score -= gapState == previousGapState ? penalties.extend() : penalties.open();
        previousGapState = gapState;
      }
    }
    return score;
  }
}
