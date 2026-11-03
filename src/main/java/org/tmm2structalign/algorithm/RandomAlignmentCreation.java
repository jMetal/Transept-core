package org.tmm2structalign.algorithm;

import java.util.ArrayList;
import java.util.List;
import org.tmm2structalign.alignment.GapIntervals;
import org.tmm2structalign.alignment.MsaSolution;
import org.tmm2structalign.model.ProteinSequence;
import org.tmm2structalign.operator.RandomSource;
import org.tmm2structalign.problem.StructuralMsaProblem;
import org.uma.jmetal.component.catalogue.common.solutionscreation.SolutionsCreation;

/**
 * Creates alignments with gaps placed uniformly at random, as a weak baseline (random search).
 *
 * <p>Every alignment has {@code ceil(longest * (1 + lengthIncrease))} columns, where {@code
 * longest} is the length of the longest protein; each row receives the gaps it needs at columns
 * chosen uniformly at random. Gap-only columns may appear; the problem removes them on
 * evaluation.
 */
public final class RandomAlignmentCreation implements SolutionsCreation<MsaSolution> {

  /** Extra length, relative to the longest protein, used by default. */
  public static final double DEFAULT_LENGTH_INCREASE = 0.2;

  private final StructuralMsaProblem problem;
  private final int batchSize;
  private final int alignmentLength;
  private final RandomSource random;

  /**
   * Creates the generator.
   *
   * @param problem problem whose proteins are aligned
   * @param batchSize number of alignments created per call
   * @param lengthIncrease extra length relative to the longest protein, non-negative
   * @param random source of random numbers
   */
  public RandomAlignmentCreation(
      StructuralMsaProblem problem, int batchSize, double lengthIncrease, RandomSource random) {
    if (batchSize < 1) {
      throw new IllegalArgumentException("Batch size must be positive: " + batchSize);
    }
    if (!(lengthIncrease >= 0.0) || Double.isInfinite(lengthIncrease)) {
      throw new IllegalArgumentException(
          "Length increase must be finite and non-negative: " + lengthIncrease);
    }
    int longest =
        problem.proteins().proteins().stream().mapToInt(ProteinSequence::length).max().orElse(0);
    this.problem = problem;
    this.batchSize = batchSize;
    this.alignmentLength = (int) Math.ceil(longest * (1.0 + lengthIncrease));
    this.random = random;
  }

  /** Returns the number of columns of the created alignments. */
  public int alignmentLength() {
    return alignmentLength;
  }

  @Override
  public List<MsaSolution> create() {
    List<MsaSolution> solutions = new ArrayList<>(batchSize);
    for (int index = 0; index < batchSize; index++) {
      List<GapIntervals> rows = new ArrayList<>(problem.numberOfVariables());
      for (ProteinSequence protein : problem.proteins().proteins()) {
        rows.add(randomRow(alignmentLength - protein.length()));
      }
      solutions.add(new MsaSolution(problem.proteins(), rows, problem.numberOfObjectives()));
    }
    return solutions;
  }

  private GapIntervals randomRow(int gaps) {
    boolean[] gapMask = new boolean[alignmentLength];
    int placed = 0;
    while (placed < gaps) {
      int column = random.nextInt(0, alignmentLength - 1);
      if (!gapMask[column]) {
        gapMask[column] = true;
        placed++;
      }
    }
    return GapIntervals.fromMask(gapMask);
  }
}
