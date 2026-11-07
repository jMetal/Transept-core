package org.transept.algorithm;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.logging.Logger;
import org.transept.alignment.Msa;
import org.transept.alignment.MsaSolution;
import org.transept.evolver.MsaConfigurations;
import org.transept.evolver.MsaGeneticAlgorithm;
import org.transept.evolver.MsaNsga2;
import org.transept.evolver.MsaOperators;
import org.transept.evolver.MsaRunContext;
import org.transept.evolver.MsaSpea2;
import org.transept.operator.MergeAdjacentGapRunsMutation;
import org.transept.operator.RandomChoiceMutation;
import org.transept.operator.RandomSource;
import org.transept.operator.TmSafeBlockShiftMutation;
import org.transept.operator.TmSafeGapShiftMutation;
import org.transept.operator.TmSafeSinglePointCrossover;
import org.transept.operator.TmSafeWindowRealignMutation;
import org.transept.problem.StructuralMsaProblem;
import org.transept.structure.ProteinStructures;
import org.uma.evolver.algorithm.BaseLevelAlgorithm;
import org.uma.evolver.parameter.ParameterSpace;
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
import org.uma.jmetal.operator.crossover.CrossoverOperator;
import org.uma.jmetal.operator.mutation.MutationOperator;
import org.uma.jmetal.operator.selection.impl.BinaryTournamentSelection;
import org.uma.jmetal.parallel.asynchronous.algorithm.impl.AsynchronousMultiThreadedGeneticAlgorithm;
import org.uma.jmetal.parallel.asynchronous.algorithm.impl.AsynchronousMultiThreadedNSGAII;
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

  /** Name of the asynchronous NSGA-II in the outputs. */
  public static final String ASYNC_NSGA_II_NAME = "AsyncNSGAII";

  /** Name of the asynchronous genetic algorithm in the outputs. */
  public static final String ASYNC_GENETIC_ALGORITHM_NAME = "AsyncGeneticAlgorithm";

  /**
   * Classpath directory of the weight vectors of MOEA/D for more than two objectives (files
   * {@code W<objectives>D_<population size>.dat}, taken from Evolver's resources).
   */
  public static final String MOEAD_WEIGHT_VECTORS = "/weightVectors";

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
    if (!configuration.configuration().isBlank() && !type.isConfigurable()) {
      throw new IllegalArgumentException(
          "A configuration string is only available for the configurable algorithms, not for "
              + type);
    }
    if (type.isSteadyState() && configuration.threads() > 1) {
      LOGGER.warning(
          type + " evaluates one solution per iteration; " + configuration.threads()
              + " threads will not speed it up");
    }
    return switch (type) {
      case NSGA_II ->
          configured(
              type,
              space ->
                  new MsaNsga2(
                      problem,
                      configuration.populationSize(),
                      configuration.maxEvaluations(),
                      space,
                      context()));
      case SPEA2 ->
          configured(
              type,
              space ->
                  new MsaSpea2(
                      problem,
                      configuration.populationSize(),
                      configuration.maxEvaluations(),
                      space,
                      context()));
      case GENETIC_ALGORITHM ->
          configured(
              type,
              space ->
                  new MsaGeneticAlgorithm(
                      problem,
                      configuration.populationSize(),
                      configuration.maxEvaluations(),
                      space,
                      context()));
      case ASYNC_NSGA_II ->
          asynchronous(
              ASYNC_NSGA_II_NAME,
              new AsynchronousMultiThreadedNSGAII<>(
                  configuration.threads(),
                  problem,
                  configuration.populationSize(),
                  crossover(),
                  mutation(),
                  termination()));
      case ASYNC_GENETIC_ALGORITHM -> {
        ObjectiveComparator<MsaSolution> fitness = new ObjectiveComparator<>(0);
        yield asynchronous(
            ASYNC_GENETIC_ALGORITHM_NAME,
            new AsynchronousMultiThreadedGeneticAlgorithm<>(
                configuration.threads(),
                problem,
                configuration.populationSize(),
                crossover(),
                mutation(),
                new BinaryTournamentSelection<>(fitness),
                new MuPlusLambdaReplacement<>(fitness),
                termination()));
      }
      case NSGA_III -> SearchRun.of(nsgaThree());
      case AGE_MOEA -> SearchRun.of(ageMoea(new AGEMOEAEnvironmentalSelection<>(objectives())));
      case AGE_MOEA_II ->
          SearchRun.of(ageMoea(new AGEMOEA2EnvironmentalSelection<>(objectives())));
      case SMS_EMOA -> SearchRun.of(smsEmoa());
      case MOEA_D -> SearchRun.of(moead());
      case PAES -> SearchRun.of(paes());
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

  /**
   * An algorithm built with Evolver from its parameter space. The configuration is resolved in
   * layers: the default configuration, the operator settings of the {@link
   * AlgorithmConfiguration} and, last, its configuration string.
   */
  @SuppressWarnings("unchecked")
  private SearchRun configured(
      AlgorithmType type, Function<ParameterSpace, BaseLevelAlgorithm<MsaSolution>> algorithm) {
    Map<String, String> values =
        new LinkedHashMap<>(MsaConfigurations.defaults(type.parameterSpace()));
    values.put(
        "offspringPopulationSize", Integer.toString(configuration.offspringPopulationSize()));
    values.put("crossoverProbability", Double.toString(configuration.crossoverProbability()));
    values.put("cutCandidates", Integer.toString(configuration.cutCandidates()));
    values.put("mutation", mutationName(configuration.mutationType()));
    values.put("mutationProbability", Double.toString(configuration.mutationProbability()));
    values.put("maxShift", Integer.toString(configuration.maxShift()));
    values.putAll(MsaConfigurations.parse(configuration.configuration()));
    ParameterSpace space = MsaConfigurations.parameterSpace(type.parameterSpace());
    MsaConfigurations.requireKnownNames(values, space);
    var built =
        (EvolutionaryAlgorithm<MsaSolution>)
            algorithm.apply(space).parse(MsaConfigurations.arguments(values)).build();
    return SearchRun.withConfiguration(SearchRun.of(built), MsaConfigurations.format(values));
  }

  private MsaRunContext context() {
    return new MsaRunContext(
        seeds, structures, random, configuration.threads(), configuration.maxTime());
  }

  /** An asynchronous algorithm of jMetal that starts from the seeded population. */
  private SearchRun asynchronous(
      String name, AsynchronousMultiThreadedGeneticAlgorithm<MsaSolution> algorithm) {
    algorithm.setInitialSolutionsCreation(population());
    return SearchRun.of(algorithm, name);
  }

  private static String mutationName(MutationType type) {
    return switch (type) {
      case TM_SAFE_GAP_SHIFT -> MsaOperators.TM_SAFE_GAP_SHIFT;
      case MERGE_GAP_RUNS -> MsaOperators.MERGE_GAP_RUNS;
      case RANDOM_CHOICE -> MsaOperators.RANDOM_CHOICE;
      case TM_SAFE_BLOCK_SHIFT -> MsaOperators.TM_SAFE_BLOCK_SHIFT;
      case TM_SAFE_WINDOW_REALIGN -> MsaOperators.TM_SAFE_WINDOW_REALIGN;
    };
  }

  private void checkObjectives(AlgorithmType type) {
    if (objectives() < type.minObjectives() || objectives() > type.maxObjectives()) {
      throw new IllegalArgumentException(
          type + " accepts " + type.minObjectives() + " to " + type.maxObjectives()
              + " objectives, got " + objectives()
              + (type.isSingleObjective()
                  ? "; use --objectives LDDT"
                  : "; use GENETIC_ALGORITHM for one objective"));
    }
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
    if (objectives() > 2) {
      String file = MOEAD_WEIGHT_VECTORS + "/W" + objectives() + "D_" + size + ".dat";
      if (AlgorithmFactory.class.getResource(file) == null) {
        throw new IllegalArgumentException(
            "MOEA_D with " + objectives() + " objectives needs a population size with weight"
                + " vectors (" + file + " is missing); with 3 objectives: 91, 100, 136, 300, 500,"
                + " 600, 800, 1000 or 1200");
      }
    }
    return new MOEADBuilder<>(
            problem,
            size,
            crossover(),
            mutation(),
            MOEAD_WEIGHT_VECTORS,
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
        problem,
        seeds,
        configuration.populationSize(),
        random,
        mutation(1.0),
        configuration.cutCandidates());
  }

  /** MOEA/D population: seeds spread over equally spaced subproblems. */
  private SolutionsCreation<MsaSolution> moeadPopulation() {
    SeededPopulationCreation creation = population();
    int seedCount = creation.distinctSeeds().size();
    return () -> SeededPopulationCreation.spreadSeeds(creation.create(), seedCount);
  }

  private CrossoverOperator<MsaSolution> crossover() {
    return new TmSafeSinglePointCrossover(
        configuration.crossoverProbability(), configuration.cutCandidates(), random);
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
              probability,
              structures,
              TmSafeGapShiftMutation.DEFAULT_MAX_TRIES,
              configuration.maxShift(),
              random);
      case MERGE_GAP_RUNS -> new MergeAdjacentGapRunsMutation(probability, random);
      case TM_SAFE_BLOCK_SHIFT ->
          new TmSafeBlockShiftMutation(
              probability, TmSafeBlockShiftMutation.DEFAULT_MAX_SHIFT, random);
      case TM_SAFE_WINDOW_REALIGN ->
          new TmSafeWindowRealignMutation(
              probability, TmSafeWindowRealignMutation.DEFAULT_WINDOW, random);
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
