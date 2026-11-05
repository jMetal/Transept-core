package org.tmm2structalign.algorithm;

import java.util.List;
import java.util.logging.Logger;
import org.tmm2structalign.alignment.Msa;
import org.tmm2structalign.alignment.MsaSolution;
import org.tmm2structalign.operator.MergeAdjacentGapRunsMutation;
import org.tmm2structalign.operator.RandomChoiceMutation;
import org.tmm2structalign.operator.RandomSource;
import org.tmm2structalign.operator.TmSafeGapShiftMutation;
import org.tmm2structalign.operator.TmSafeSinglePointCrossover;
import org.tmm2structalign.problem.StructuralMsaProblem;
import org.tmm2structalign.structure.ProteinStructures;
import org.uma.jmetal.component.algorithm.EvolutionaryAlgorithm;
import org.uma.jmetal.component.algorithm.RandomSearchAlgorithm;
import org.uma.jmetal.component.algorithm.multiobjective.AGEMOEABuilder;
import org.uma.jmetal.component.algorithm.multiobjective.MOEADBuilder;
import org.uma.jmetal.component.algorithm.multiobjective.NSGAIIIBuilder;
import org.uma.jmetal.component.algorithm.multiobjective.PAESBuilder;
import org.uma.jmetal.component.algorithm.multiobjective.SMSEMOABuilder;
import org.uma.jmetal.component.catalogue.common.evaluation.Evaluation;
import org.uma.jmetal.component.catalogue.common.evaluation.impl.MultiThreadedEvaluation;
import org.uma.jmetal.component.catalogue.common.evaluation.impl.SequentialEvaluation;
import org.uma.jmetal.component.catalogue.common.solutionscreation.SolutionsCreation;
import org.uma.jmetal.component.catalogue.common.termination.Termination;
import org.uma.jmetal.component.catalogue.common.termination.impl.TerminationByComputingTime;
import org.uma.jmetal.component.catalogue.common.termination.impl.TerminationByEvaluations;
import org.uma.jmetal.component.catalogue.ea.replacement.impl.MuPlusLambdaReplacement;
import org.uma.jmetal.component.catalogue.ea.replacement.impl.agemoea.AGEMOEA2EnvironmentalSelection;
import org.uma.jmetal.component.catalogue.ea.replacement.impl.agemoea.AGEMOEAEnvironmentalSelection;
import org.uma.jmetal.component.catalogue.ea.selection.impl.NaryTournamentSelection;
import org.uma.jmetal.component.catalogue.ea.variation.impl.CrossoverAndMutationVariation;
import org.uma.jmetal.operator.crossover.CrossoverOperator;
import org.uma.jmetal.operator.mutation.MutationOperator;
import org.uma.jmetal.util.archive.impl.CrowdingDistanceArchive;
import org.uma.jmetal.util.comparator.ObjectiveComparator;
import org.uma.jmetal.util.sequencegenerator.impl.CyclicIntegerSequence;

/**
 * Builds every {@link AlgorithmType} with the project's seeded initial population and operators.
 *
 * <p>All population-based algorithms start from {@link SeededPopulationCreation}; PAES starts
 * from the first distinct seed; random search does not use seeds ({@code RANDOM_SEARCH}) or only
 * perturbs them ({@code RANDOM_SEEDED}). Every algorithm stops after the configured number of
 * evaluations.
 */
public final class AlgorithmFactory {

  /** Name reported by jMetal for the genetic algorithm. */
  public static final String GENETIC_ALGORITHM_NAME = "GeneticAlgorithm";

  /** Neighbourhood size of MOEA/D (capped by the population size). */
  public static final int MOEAD_NEIGHBORHOOD_SIZE = 20;

  private static final Logger LOGGER = Logger.getLogger(AlgorithmFactory.class.getName());

  private final StructuralMsaProblem problem;
  private final List<Msa> seeds;
  private final AlgorithmConfiguration configuration;
  private final ProteinStructures structures;
  private final RandomSource random;

  /**
   * Creates the factory.
   *
   * @param problem problem to solve
   * @param seeds seed alignments in order of preference
   * @param configuration algorithm parameters
   * @param structures structures used by the mutation
   * @param random source of random numbers of the operators and the initial population
   */
  public AlgorithmFactory(
      StructuralMsaProblem problem,
      List<Msa> seeds,
      AlgorithmConfiguration configuration,
      ProteinStructures structures,
      RandomSource random) {
    this.problem = problem;
    this.seeds = List.copyOf(seeds);
    this.configuration = configuration;
    this.structures = structures;
    this.random = random;
  }

  /**
   * Builds an algorithm.
   *
   * @param type algorithm to build
   * @return the run, ready to start
   * @throws IllegalArgumentException if the configuration does not suit the algorithm, or the
   *     number of objectives does not suit it: the genetic algorithm needs exactly one, the
   *     multiobjective algorithms at least two (the random searches accept any number)
   */
  public SearchRun create(AlgorithmType type) {
    checkObjectives(type);
    if (configuration.externalArchive() && type != AlgorithmType.NSGA_II) {
      throw new IllegalArgumentException("The external archive is only available with NSGA_II");
    }
    if (type.isSteadyState() && configuration.threads() > 1) {
      LOGGER.warning(
          type + " evaluates one solution per iteration; " + configuration.threads()
              + " threads will not speed it up");
    }
    return switch (type) {
      case NSGA_II -> nsgaTwo();
      case NSGA_III -> SearchRun.of(nsgaThree());
      case AGE_MOEA -> SearchRun.of(ageMoea(new AGEMOEAEnvironmentalSelection<>(objectives())));
      case AGE_MOEA_II ->
          SearchRun.of(ageMoea(new AGEMOEA2EnvironmentalSelection<>(objectives())));
      case SMS_EMOA -> SearchRun.of(smsEmoa());
      case MOEA_D -> SearchRun.of(moead());
      case PAES -> SearchRun.of(paes());
      case GENETIC_ALGORITHM -> SearchRun.of(geneticAlgorithm());
      case RANDOM_SEARCH ->
          SearchRun.of(
              new RandomSearchAlgorithm<>(
                  "RandomSearch",
                  new RandomAlignmentCreation(
                      problem,
                      configuration.populationSize(),
                      RandomAlignmentCreation.DEFAULT_LENGTH_INCREASE,
                      random),
                  evaluation(),
                  termination()));
      case RANDOM_SEEDED ->
          SearchRun.of(
              new RandomSearchAlgorithm<>(
                  "RandomSeededSearch",
                  new PerturbedSeedCreation(
                      problem,
                      population().distinctSeeds(),
                      configuration.populationSize(),
                      PerturbedSeedCreation.DEFAULT_MUTATIONS,
                      mutation(1.0),
                      random),
                  evaluation(),
                  termination()));
    };
  }

  /** NSGA-II whose result is the final population or, if configured, a bounded archive. */
  private SearchRun nsgaTwo() {
    if (!configuration.externalArchive()) {
      return SearchRun.of(
          TmStructAlign.create(problem, population(), configuration, crossover(), mutation()));
    }
    var archive = new CrowdingDistanceArchive<MsaSolution>(configuration.populationSize());
    var algorithm =
        TmStructAlign.create(
            problem,
            population(),
            configuration,
            crossover(),
            mutation(),
            new ArchivingEvaluation(evaluation(), archive));
    return SearchRun.of(algorithm, archive);
  }

  private void checkObjectives(AlgorithmType type) {
    boolean randomSearch =
        type == AlgorithmType.RANDOM_SEARCH || type == AlgorithmType.RANDOM_SEEDED;
    if (type == AlgorithmType.GENETIC_ALGORITHM && objectives() != 1) {
      throw new IllegalArgumentException(
          "GENETIC_ALGORITHM optimises a single objective, got " + objectives()
              + "; use --objectives LDDT");
    }
    if (type != AlgorithmType.GENETIC_ALGORITHM && !randomSearch && objectives() < 2) {
      throw new IllegalArgumentException(
          type + " needs at least two objectives; use GENETIC_ALGORITHM for one");
    }
  }

  /**
   * Generational genetic algorithm: binary tournament on the single objective, the same
   * variation as NSGA-II and (mu + lambda) replacement, which keeps the best solutions.
   */
  private EvolutionaryAlgorithm<MsaSolution> geneticAlgorithm() {
    var variation =
        new CrossoverAndMutationVariation<>(
            configuration.offspringPopulationSize(), crossover(), mutation());
    ObjectiveComparator<MsaSolution> fitness = new ObjectiveComparator<>(0);
    return new EvolutionaryAlgorithm<>(
        GENETIC_ALGORITHM_NAME,
        population(),
        evaluation(),
        termination(),
        new NaryTournamentSelection<>(2, variation.matingPoolSize(), fitness),
        variation,
        new MuPlusLambdaReplacement<>(fitness));
  }

  private EvolutionaryAlgorithm<MsaSolution> nsgaThree() {
    int size = configuration.populationSize();
    if (size % 4 != 0 || size < 4) {
      throw new IllegalArgumentException(
          "NSGA-III needs a population size that is a multiple of 4, got " + size);
    }
    return new NSGAIIIBuilder<>(problem, size - 1, crossover(), mutation())
        .setCreateInitialPopulation(population())
        .setTermination(termination())
        .setEvaluation(evaluation())
        .build();
  }

  private EvolutionaryAlgorithm<MsaSolution> ageMoea(
      AGEMOEAEnvironmentalSelection<MsaSolution> environmentalSelection) {
    return new AGEMOEABuilder<>(
            problem,
            configuration.populationSize(),
            configuration.offspringPopulationSize(),
            crossover(),
            mutation(),
            environmentalSelection)
        .setCreateInitialPopulation(population())
        .setTermination(termination())
        .setEvaluation(evaluation())
        .build();
  }

  private EvolutionaryAlgorithm<MsaSolution> smsEmoa() {
    return new SMSEMOABuilder<>(problem, configuration.populationSize(), crossover(), mutation())
        .setCreateInitialPopulation(population())
        .setTermination(termination())
        .setEvaluation(evaluation())
        .build();
  }

  private EvolutionaryAlgorithm<MsaSolution> moead() {
    int size = configuration.populationSize();
    return new MOEADBuilder<>(
            problem,
            size,
            crossover(),
            mutation(),
            "",
            new CyclicIntegerSequence(size),
            true)
        .setNeighborhoodSize(Math.min(MOEAD_NEIGHBORHOOD_SIZE, size))
        .setCreateInitialPopulation(moeadPopulation())
        .setTermination(termination())
        .setEvaluation(evaluation())
        .build();
  }

  private EvolutionaryAlgorithm<MsaSolution> paes() {
    int archiveSize = configuration.populationSize();
    MsaSolution initialSolution = problem.createSolution(population().distinctSeeds().get(0));
    return new PAESBuilder<>(
            problem,
            archiveSize,
            mutation(),
            new CrowdingDistanceArchive<>(archiveSize),
            initialSolution)
        .setTermination(termination())
        .setEvaluation(new SequentialEvaluation<>(problem))
        .build();
  }

  private SeededPopulationCreation population() {
    return new SeededPopulationCreation(
        problem, seeds, configuration.populationSize(), random, mutation(1.0));
  }

  /** MOEA/D population: seeds spread over equally spaced subproblems. */
  private SolutionsCreation<MsaSolution> moeadPopulation() {
    SeededPopulationCreation creation = population();
    int seedCount = creation.distinctSeeds().size();
    return () -> SeededPopulationCreation.spreadSeeds(creation.create(), seedCount);
  }

  private CrossoverOperator<MsaSolution> crossover() {
    return new TmSafeSinglePointCrossover(configuration.crossoverProbability(), random);
  }

  private MutationOperator<MsaSolution> mutation() {
    return mutation(configuration.mutationProbability());
  }

  private MutationOperator<MsaSolution> mutation(double probability) {
    return mutation(configuration.mutationType(), probability);
  }

  private MutationOperator<MsaSolution> mutation(MutationType type, double probability) {
    return switch (type) {
      case TM_SAFE_GAP_SHIFT ->
          new TmSafeGapShiftMutation(
              probability, structures, TmSafeGapShiftMutation.DEFAULT_MAX_TRIES, random);
      case MERGE_GAP_RUNS -> new MergeAdjacentGapRunsMutation(probability, random);
      case RANDOM_CHOICE ->
          new RandomChoiceMutation(
              probability,
              List.of(
                  mutation(MutationType.TM_SAFE_GAP_SHIFT, 1.0),
                  mutation(MutationType.MERGE_GAP_RUNS, 1.0)),
              1,
              random);
    };
  }

  private Termination termination() {
    return configuration.maxTime() == null
        ? new TerminationByEvaluations(configuration.maxEvaluations())
        : new TerminationByComputingTime(configuration.maxTime().toMillis());
  }

  private Evaluation<MsaSolution> evaluation() {
    return configuration.threads() > 1
        ? new MultiThreadedEvaluation<>(configuration.threads(), problem)
        : new SequentialEvaluation<>(problem);
  }

  private int objectives() {
    return problem.numberOfObjectives();
  }
}
