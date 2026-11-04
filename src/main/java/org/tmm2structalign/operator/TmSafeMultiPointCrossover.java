package org.tmm2structalign.operator;

import java.util.List;
import org.tmm2structalign.alignment.MsaSolution;
import org.uma.jmetal.operator.crossover.CrossoverOperator;

/**
 * Crossover with several cuts, none of which splits a transmembrane helix.
 *
 * <p>Generalises {@link TmSafeSinglePointCrossover}: the child is built by repeated single-point
 * recombinations that alternate the donor parent. With two points, the child takes the
 * beginning of the first parent, a middle segment of the second and the end of the first. Each
 * cut is chosen uniformly among the valid cuts of the partial child that lie after the previous
 * cut, so the columns already copied are kept. If no further valid cut exists the child has
 * fewer cuts; if there is none at all it is a copy of the first parent.
 *
 * <p>Every crossover produces two children, (first, second) and (second, first), each with its
 * own probability draw.
 */
public final class TmSafeMultiPointCrossover implements CrossoverOperator<MsaSolution> {

  private final double probability;
  private final int points;
  private final RandomSource random;
  private final TmSafeSinglePointCrossover singlePoint;

  /**
   * Creates the operator using the global jMetal generator.
   *
   * @param probability probability of applying the crossover to each child
   * @param points number of cuts, at least 1
   */
  public TmSafeMultiPointCrossover(double probability, int points) {
    this(probability, points, RandomSource.jmetal());
  }

  /**
   * Creates the operator.
   *
   * @param probability probability of applying the crossover to each child
   * @param points number of cuts, at least 1
   * @param random source of random numbers
   */
  public TmSafeMultiPointCrossover(double probability, int points, RandomSource random) {
    if (points < 1) {
      throw new IllegalArgumentException("The number of points must be positive: " + points);
    }
    this.singlePoint = new TmSafeSinglePointCrossover(probability, random);
    this.probability = probability;
    this.points = points;
    this.random = random;
  }

  @Override
  public List<MsaSolution> execute(List<MsaSolution> parents) {
    if (parents == null || parents.size() != 2) {
      throw new IllegalArgumentException("The crossover needs exactly two parents");
    }
    MsaSolution first = parents.get(0);
    MsaSolution second = parents.get(1);
    if (first.proteins() != second.proteins()) {
      throw new IllegalArgumentException("Parents must align the same protein set");
    }
    return List.of(offspring(first, second), offspring(second, first));
  }

  @Override
  public double crossoverProbability() {
    return probability;
  }

  @Override
  public int numberOfRequiredParents() {
    return 2;
  }

  @Override
  public int numberOfGeneratedChildren() {
    return 2;
  }

  /** Returns the number of cuts. */
  public int points() {
    return points;
  }

  /**
   * Builds the child of a list of cuts, applied in order.
   *
   * @param first parent providing the first segment
   * @param second parent providing the second segment
   * @param cuts increasing cut columns, each valid in the partial child built so far
   * @return the child, not evaluated
   */
  public MsaSolution recombine(MsaSolution first, MsaSolution second, List<Integer> cuts) {
    MsaSolution child = first.copy();
    for (int index = 0; index < cuts.size(); index++) {
      child = singlePoint.recombine(child, index % 2 == 0 ? second : first, cuts.get(index));
    }
    return child;
  }

  private MsaSolution offspring(MsaSolution first, MsaSolution second) {
    MsaSolution child = first.copy();
    if (random.nextDouble() < probability) {
      int previousCut = 0;
      for (int index = 0; index < points; index++) {
        final int after = previousCut;
        List<Integer> cuts =
            singlePoint.compatibleCuts(child).stream().filter(cut -> cut > after).toList();
        if (cuts.isEmpty()) {
          break;
        }
        int cut = cuts.get(random.nextInt(0, cuts.size() - 1));
        child = singlePoint.recombine(child, index % 2 == 0 ? second : first, cut);
        previousCut = cut;
      }
    }
    return child;
  }
}
