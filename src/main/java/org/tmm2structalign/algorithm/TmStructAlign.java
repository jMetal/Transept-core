package org.tmm2structalign.algorithm;

import java.util.Comparator;
import java.util.List;
import org.tmm2structalign.alignment.MsaSolution;
import org.tmm2structalign.operator.RandomSource;
import org.tmm2structalign.operator.TmSafeSinglePointCrossover;
import org.tmm2structalign.problem.StructuralMsaProblem;
import org.uma.jmetal.component.algorithm.EvolutionaryAlgorithm;
import org.uma.jmetal.component.catalogue.common.evaluation.Evaluation;
import org.uma.jmetal.component.catalogue.common.evaluation.impl.MultiThreadedEvaluation;
import org.uma.jmetal.component.catalogue.common.evaluation.impl.SequentialEvaluation;
import org.uma.jmetal.component.catalogue.common.termination.impl.TerminationByEvaluations;
import org.uma.jmetal.component.catalogue.ea.replacement.impl.RankingAndDensityEstimatorReplacement;
import org.uma.jmetal.component.catalogue.ea.selection.impl.NaryTournamentSelection;
import org.uma.jmetal.component.catalogue.ea.variation.impl.CrossoverAndMutationVariation;
import org.uma.jmetal.operator.crossover.CrossoverOperator;
import org.uma.jmetal.operator.mutation.MutationOperator;
import org.uma.jmetal.util.comparator.MultiComparator;
import org.uma.jmetal.util.densityestimator.impl.CrowdingDistanceDensityEstimator;
import org.uma.jmetal.util.ranking.impl.FastNonDominatedSortRanking;

/**
 * Factory of the TM-M2StructAlign algorithm: an NSGA-II built from jMetal components, with
 * binary tournament selection (rank, then crowding distance), the given variation operators and
 * ranking-and-crowding replacement.
 */
public final class TmStructAlign {

  /** Algorithm name reported by jMetal. */
  public static final String NAME = "TM-M2StructAlign";

  private static final int TOURNAMENT_SIZE = 2;

  private TmStructAlign() {}

  /**
   * Creates the algorithm.
   *
   * @param problem problem to solve
   * @param population initial population
   * @param configuration algorithm parameters
   * @param crossover crossover operator
   * @param mutation mutation operator
   * @return the algorithm, ready to run
   */
  public static EvolutionaryAlgorithm<MsaSolution> create(
      StructuralMsaProblem problem,
      SeededPopulationCreation population,
      AlgorithmConfiguration configuration,
      CrossoverOperator<MsaSolution> crossover,
      MutationOperator<MsaSolution> mutation) {
    Evaluation<MsaSolution> evaluation =
        configuration.threads() > 1
            ? new MultiThreadedEvaluation<>(configuration.threads(), problem)
            : new SequentialEvaluation<>(problem);
    return create(problem, population, configuration, crossover, mutation, evaluation);
  }

  /**
   * Creates the algorithm with a given evaluation, e.g. one that also fills an external archive.
   *
   * @param problem problem to solve
   * @param population initial population
   * @param configuration algorithm parameters
   * @param crossover crossover operator
   * @param mutation mutation operator
   * @param evaluation evaluation of the solutions
   * @return the algorithm, ready to run
   */
  public static EvolutionaryAlgorithm<MsaSolution> create(
      StructuralMsaProblem problem,
      SeededPopulationCreation population,
      AlgorithmConfiguration configuration,
      CrossoverOperator<MsaSolution> crossover,
      MutationOperator<MsaSolution> mutation,
      Evaluation<MsaSolution> evaluation) {
    var variation =
        new CrossoverAndMutationVariation<>(
            configuration.offspringPopulationSize(), crossover, mutation);
    var ranking = new FastNonDominatedSortRanking<MsaSolution>();
    var density = new CrowdingDistanceDensityEstimator<MsaSolution>();
    var selection =
        new NaryTournamentSelection<>(
            TOURNAMENT_SIZE,
            variation.matingPoolSize(),
            new MultiComparator<MsaSolution>(
                List.of(
                    Comparator.comparing(ranking::getRank),
                    Comparator.comparing(density::value).reversed())));
    var replacement =
        new RankingAndDensityEstimatorReplacement<>(
            ranking, density, RankingAndDensityEstimatorReplacement.RemovalPolicy.ONE_SHOT);
    return new EvolutionaryAlgorithm<>(
        NAME,
        population,
        evaluation,
        new TerminationByEvaluations(configuration.maxEvaluations()),
        selection,
        variation,
        replacement);
  }

  /**
   * Creates the algorithm with the published operators.
   *
   * @param problem problem to solve
   * @param population initial population
   * @param configuration algorithm parameters
   * @param mutation mutation operator (needs the problem's structures)
   * @param random source of random numbers for the crossover
   * @return the algorithm, ready to run
   */
  public static EvolutionaryAlgorithm<MsaSolution> createDefault(
      StructuralMsaProblem problem,
      SeededPopulationCreation population,
      AlgorithmConfiguration configuration,
      MutationOperator<MsaSolution> mutation,
      RandomSource random) {
    return create(
        problem,
        population,
        configuration,
        new TmSafeSinglePointCrossover(configuration.crossoverProbability(), random),
        mutation);
  }
}
