package org.transept.algorithm;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assumptions.assumeTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Duration;
import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.EnumSource;
import org.transept.alignment.Msa;
import org.transept.alignment.MsaSolution;
import org.transept.evolver.MsaGeneticAlgorithm;
import org.transept.model.ProteinSet;
import org.transept.objective.AlignedPairCoverage;
import org.transept.objective.StructuralLddt;
import org.transept.objective.TopologyAwareSumOfPairs;
import org.transept.problem.StructuralMsaProblem;
import org.transept.structure.ProteinStructures;
import org.transept.testing.Alignments;
import org.transept.testing.SeededRandom;

@DisplayName("Unit tests for classes AlgorithmFactory, SearchRun and AlgorithmType")
class AlgorithmFactoryTest {

  private static final int POPULATION_SIZE = 8;
  private static final int EVALUATIONS = 40;

  private ProteinSet proteins;
  private StructuralMsaProblem problem;
  private ProteinStructures structures;
  private List<Msa> seeds;

  @BeforeEach
  void setUp() {
    proteins =
        Alignments.proteins("ACDEFGHIK:OOMMMMOOO", "ACDEGHIK:OOMMMOOO", "ACEFGHIK:OOMMMOOO");
    structures = Alignments.fullyStructured(proteins);
    problem =
        new StructuralMsaProblem(
            "test",
            proteins,
            List.of(
                TopologyAwareSumOfPairs.withDefaultParameters(),
                new StructuralLddt(structures, 15f)));
    seeds =
        List.of(
            Msa.fromAlignedRows(proteins, List.of("ACDEFGHIK-", "ACDE-GHIK-", "AC-EFGHIK-")),
            Msa.fromAlignedRows(proteins, List.of("ACDEFGHIK", "ACD-EGHIK", "A-CEFGHIK")),
            Msa.fromAlignedRows(proteins, List.of("-ACDEFGHIK", "-ACDEGHIK-", "--ACEFGHIK")));
  }

  private AlgorithmFactory factory(int populationSize, double mutationProbability) {
    return factory(problem, populationSize, mutationProbability);
  }

  /** Factory whose problem suits the algorithm: one objective for the genetic algorithm. */
  private AlgorithmFactory factory(
      AlgorithmType type, int populationSize, double mutationProbability) {
    StructuralMsaProblem suited =
        type.isSingleObjective()
            ? new StructuralMsaProblem(
                "test", proteins, List.of(new StructuralLddt(structures, 15f)))
            : problem;
    return factory(suited, populationSize, mutationProbability);
  }

  private AlgorithmFactory factory(
      StructuralMsaProblem chosen, int populationSize, double mutationProbability) {
    return new AlgorithmFactory(
        chosen,
        seeds,
        new AlgorithmConfiguration(EVALUATIONS, populationSize, populationSize, 0.8,
            mutationProbability, 1),
        structures,
        new SeededRandom(5));
  }

  @Nested
  @DisplayName("When building and running every algorithm")
  class Building {

    @ParameterizedTest
    @CsvSource({
      "NSGA_II, NSGAII",
      "NSGA_III, NSGAIII",
      "AGE_MOEA, AGE-MOEA",
      "AGE_MOEA_II, AGE-MOEA-II",
      "SMS_EMOA, SMS-EMOA",
      "MOEA_D, MOEAD",
      "PAES, PAES",
      "RANDOM_SEARCH, RandomSearch",
      "RANDOM_SEEDED, RandomSeededSearch",
      "GENETIC_ALGORITHM, GeneticAlgorithm",
      "SPEA2, RDEMOEA",
      "ASYNC_NSGA_II, AsyncNSGAII",
      "ASYNC_GENETIC_ALGORITHM, AsyncGeneticAlgorithm"
    })
    @DisplayName("given an algorithm type, when built, then jMetal reports the expected name")
    void givenAlgorithmType_whenBuilt_thenJmetalReportsExpectedName(
        AlgorithmType type, String name) {
      // Act
      SearchRun run = factory(type, POPULATION_SIZE, 0.2).create(type);

      // Assert
      assertEquals(name, run.name());
    }

    @ParameterizedTest
    @EnumSource(AlgorithmType.class)
    @DisplayName("given an algorithm type, when run, then the budget is used and progress reported")
    void givenAlgorithmType_whenRun_thenBudgetIsUsedAndProgressReported(AlgorithmType type) {
      // Arrange
      SearchRun run = factory(type, POPULATION_SIZE, 1.0).create(type);
      List<Integer> progress = new ArrayList<>();
      run.onProgress(progress::add);

      // Act
      run.run();

      // Assert
      assertTrue(run.evaluations() >= EVALUATIONS, type + ": " + run.evaluations());
      assertFalse(run.front().isEmpty());
      assertFalse(progress.isEmpty());
      assertTrue(run.runtimeMillis() >= 0);
      for (MsaSolution solution : run.front()) {
        assertEquals(3, solution.toMsa().size());
      }
    }

    @Test
    @DisplayName("given PAES, when run, then its archive never exceeds the population size")
    void givenPaes_whenRun_thenArchiveNeverExceedsPopulationSize() {
      // Arrange
      SearchRun run = factory(POPULATION_SIZE, 1.0).create(AlgorithmType.PAES);

      // Act
      run.run();

      // Assert
      assertTrue(run.front().size() <= POPULATION_SIZE);
    }
  }

  @Nested
  @DisplayName("When choosing the operators")
  class Operators {

    @ParameterizedTest
    @EnumSource(MutationType.class)
    @DisplayName("given a mutation, when running NSGA-II, then the front is valid")
    void givenMutation_whenRunningNsgaTwo_thenFrontIsValid(MutationType mutation) {
      // Arrange
      AlgorithmFactory factory =
          new AlgorithmFactory(
              problem,
              seeds,
              new AlgorithmConfiguration(
                  EVALUATIONS, POPULATION_SIZE, POPULATION_SIZE, 0.9, 1.0, 1, mutation),
              structures,
              new SeededRandom(7));
      SearchRun run = factory.create(AlgorithmType.NSGA_II);

      // Act
      run.run();

      // Assert
      assertFalse(run.front().isEmpty(), mutation.toString());
      for (MsaSolution solution : run.front()) {
        assertEquals(3, solution.toMsa().size());
      }
    }

    @ParameterizedTest
    @EnumSource(
        value = AlgorithmType.class,
        names = {"NSGA_II", "SPEA2", "SMS_EMOA", "PAES", "RANDOM_SEEDED", "ASYNC_NSGA_II"})
    @DisplayName("given a stop request, when running, then it ends early with a valid front")
    void givenStopRequest_whenRunning_thenItEndsEarlyWithValidFront(
        AlgorithmType type, @TempDir Path directory) throws IOException {
      // Arrange: the stop is requested before the run starts
      Path file = Files.writeString(directory.resolve(StopRequest.FILE_NAME), "");
      StopRequest stop = StopRequest.byFile(file);
      SearchRun run = factory(POPULATION_SIZE, 0.2).stoppableBy(stop).create(type);

      // Act
      run.run();

      // Assert
      assertTrue(run.evaluations() < EVALUATIONS, type + ": " + run.evaluations());
      assertFalse(run.front().isEmpty(), type + " after " + run.evaluations() + " evaluations");
      assertTrue(stop.wasRequested());
    }

    @ParameterizedTest
    @EnumSource(
        value = AlgorithmType.class,
        names = {"NSGA_II", "GENETIC_ALGORITHM", "SMS_EMOA", "MOEA_D", "NSGA_III"})
    @DisplayName("given the uniform crossover, when running, then the front is valid")
    void givenUniformCrossover_whenRunning_thenFrontIsValid(AlgorithmType type) {
      // Arrange
      AlgorithmFactory factory =
          new AlgorithmFactory(
              type == AlgorithmType.GENETIC_ALGORITHM
                  ? new StructuralMsaProblem(
                      "test", proteins, List.of(new StructuralLddt(structures, 15f)))
                  : problem,
              seeds,
              new AlgorithmConfiguration(
                  EVALUATIONS, POPULATION_SIZE, POPULATION_SIZE, 1.0, 0.2, 1,
                  MutationType.TM_SAFE_GAP_SHIFT, null, "", 0, 1, CrossoverType.TM_SAFE_UNIFORM),
              structures,
              new SeededRandom(7));
      SearchRun run = factory.create(type);

      // Act
      run.run();

      // Assert
      assertFalse(run.front().isEmpty(), type.toString());
      if (type.isConfigurable()) {
        assertTrue(run.configuration().contains("--crossover tmSafeUniform"), run.configuration());
      }
    }
  }

  @Nested
  @DisplayName("When the configuration does not suit the algorithm")
  class InvalidConfiguration {

    @Test
    @DisplayName("given a population not multiple of 4, when building NSGA-III, then it fails")
    void givenPopulationNotMultipleOfFour_whenBuildingNsgaThree_thenItFails() {
      // Arrange
      AlgorithmFactory factory = factory(6, 0.2);

      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> factory.create(AlgorithmType.NSGA_III));
    }

    @Test
    @DisplayName("given two threads and a steady-state algorithm, when run, then it still works")
    void givenTwoThreadsAndSteadyStateAlgorithm_whenBuilding_thenItStillWorks() {
      // Arrange
      AlgorithmFactory factory =
          new AlgorithmFactory(
              problem,
              seeds,
              new AlgorithmConfiguration(
                  EVALUATIONS, POPULATION_SIZE, POPULATION_SIZE, 0.8, 0.2, 2),
              structures,
              new SeededRandom(5));

      // Act
      SearchRun run = factory.create(AlgorithmType.SMS_EMOA);
      run.run();

      // Assert
      assertTrue(run.evaluations() >= EVALUATIONS);
    }
  }

  @Nested
  @DisplayName("When describing algorithm types")
  class Types {

    @Test
    @DisplayName("given the algorithm types, when queried, then their properties are consistent")
    void givenAlgorithmTypes_whenQueried_thenPropertiesAreConsistent() {
      // Act & Assert
      assertTrue(AlgorithmType.SMS_EMOA.isSteadyState());
      assertTrue(AlgorithmType.MOEA_D.isSteadyState());
      assertFalse(AlgorithmType.NSGA_II.isSteadyState());
      assertFalse(AlgorithmType.PAES.usesCrossover());
      assertTrue(AlgorithmType.NSGA_III.usesCrossover());
      assertFalse(AlgorithmType.RANDOM_SEARCH.usesCrossover());
      assertTrue(AlgorithmType.GENETIC_ALGORITHM.usesCrossover());
    }
  }

  @Nested
  @DisplayName("When the number of objectives must suit the algorithm")
  class Formulations {

    private StructuralMsaProblem lddtOnly() {
      return new StructuralMsaProblem(
          "test", proteins, List.of(new StructuralLddt(structures, 15f)));
    }

    @Test
    @DisplayName("given two objectives, when building the genetic algorithm, then it fails")
    void givenTwoObjectives_whenBuildingGeneticAlgorithm_thenIllegalArgumentExceptionIsThrown() {
      // Arrange
      AlgorithmFactory factory = factory(POPULATION_SIZE, 0.2);

      // Act & Assert
      assertThrows(
          IllegalArgumentException.class, () -> factory.create(AlgorithmType.GENETIC_ALGORITHM));
    }

    @Test
    @DisplayName("given one objective, when building NSGA-II, then it fails")
    void givenOneObjective_whenBuildingNsgaTwo_thenIllegalArgumentExceptionIsThrown() {
      // Arrange
      AlgorithmFactory factory = factory(lddtOnly(), POPULATION_SIZE, 0.2);

      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> factory.create(AlgorithmType.NSGA_II));
    }

    @Test
    @DisplayName("given one objective, when running the seeded random search, then it works")
    void givenOneObjective_whenRunningSeededRandomSearch_thenItWorks() {
      // Arrange
      SearchRun run =
          factory(lddtOnly(), POPULATION_SIZE, 1.0).create(AlgorithmType.RANDOM_SEEDED);

      // Act
      run.run();

      // Assert
      assertFalse(run.front().isEmpty());
    }

    @Test
    @DisplayName("given the genetic algorithm, when run, then it keeps or improves the best seed")
    void givenGeneticAlgorithm_whenRun_thenItKeepsTheBestSeedOrImprovesIt() {
      // Arrange
      StructuralMsaProblem single = lddtOnly();
      double bestSeed = Double.POSITIVE_INFINITY;
      for (Msa seed : seeds) {
        MsaSolution solution = single.createSolution(seed);
        single.evaluate(solution);
        bestSeed = Math.min(bestSeed, solution.objectives()[0]);
      }
      SearchRun run =
          factory(single, POPULATION_SIZE, 0.2).create(AlgorithmType.GENETIC_ALGORITHM);

      // Act
      run.run();

      // Assert: objectives are stored minimised
      assertEquals(MsaGeneticAlgorithm.NAME, run.name());
      double best =
          run.front().stream()
              .mapToDouble(solution -> solution.objectives()[0])
              .min()
              .orElseThrow();
      assertTrue(best <= bestSeed, best + " > " + bestSeed);
    }
  }

  @Nested
  @DisplayName("When configuring NSGA-II from its parameter space")
  class Configuration {

    private AlgorithmFactory configured(String text) {
      return new AlgorithmFactory(
          problem,
          seeds,
          new AlgorithmConfiguration(
              EVALUATIONS, POPULATION_SIZE, POPULATION_SIZE, 0.8, 0.2, 1,
              MutationType.TM_SAFE_GAP_SHIFT, null, text),
          structures,
          new SeededRandom(5));
    }

    @Test
    @DisplayName("given the default options, when built, then the default configuration is used")
    void givenDefaultOptions_whenBuilt_thenTheDefaultConfigurationIsUsed() {
      // Act
      SearchRun run = configured("").create(AlgorithmType.NSGA_II);

      // Assert: the offspring size follows the population size
      assertEquals(
          "--algorithmResult population --createInitialSolutions seededDistinct"
              + " --offspringPopulationSize " + POPULATION_SIZE
              + " --variation crossoverAndMutationVariation --crossover tmSafeSinglePoint"
              + " --crossoverProbability 0.8 --cutCandidates 1"
              + " --mutation tmSafeGapShift"
              + " --mutationProbability 0.2 --maxShift 0 --blockShift 1 --realignWindow 50"
              + " --selection tournament"
              + " --selectionTournamentSize 2",
          run.configuration());
    }

    @Test
    @DisplayName("given a configuration string, when built, then it overrides the defaults")
    void givenConfigurationString_whenBuilt_thenItOverridesTheDefaults() {
      // Act
      SearchRun run =
          configured("--selection random --mutation mergeGapRuns").create(AlgorithmType.NSGA_II);
      run.run();

      // Assert
      assertTrue(run.configuration().contains("--selection random"), run.configuration());
      assertTrue(run.configuration().contains("--mutation mergeGapRuns"), run.configuration());
      assertTrue(run.configuration().contains("--crossoverProbability 0.8"), run.configuration());
      assertFalse(run.front().isEmpty());
    }

    @Test
    @DisplayName("given invalid configurations, when built, then an exception is thrown")
    void givenInvalidConfigurations_whenBuilt_thenExceptionIsThrown() {
      // Act & Assert
      assertThrows(
          IllegalArgumentException.class,
          () -> configured("--selectoin random").create(AlgorithmType.NSGA_II));
      assertThrows(
          RuntimeException.class,
          () -> configured("--selection roulette").create(AlgorithmType.NSGA_II));
      assertThrows(
          IllegalArgumentException.class,
          () -> configured("--selection random").create(AlgorithmType.SMS_EMOA));
    }
  }

  @Nested
  @DisplayName("When optimising three objectives")
  class ThreeObjectives {

    @ParameterizedTest
    @EnumSource(AlgorithmType.class)
    @DisplayName("given an algorithm that accepts three objectives, when run, then it completes")
    void givenAlgorithmAcceptingThreeObjectives_whenRun_thenItCompletes(AlgorithmType type) {
      // Arrange: MOEA/D needs a population size with weight vectors (W3D_100.dat)
      assumeTrue(type.maxObjectives() == 3, type + " accepts fewer objectives");
      StructuralMsaProblem three =
          new StructuralMsaProblem(
              "test",
              proteins,
              List.of(
                  TopologyAwareSumOfPairs.withDefaultParameters(),
                  new StructuralLddt(structures, 15f),
                  new AlignedPairCoverage(proteins)));
      SearchRun run = factory(three, 100, 1.0).create(type);

      // Act
      run.run();

      // Assert
      assertFalse(run.front().isEmpty());
      run.front().forEach(solution -> assertEquals(3, solution.objectives().length));
    }

    @Test
    @DisplayName("given MOEA/D, three objectives and no weights for the size, then it fails")
    void givenMoeadThreeObjectivesAndNoWeights_whenBuilding_thenItFails() {
      // Arrange
      StructuralMsaProblem three =
          new StructuralMsaProblem(
              "test",
              proteins,
              List.of(
                  TopologyAwareSumOfPairs.withDefaultParameters(),
                  new StructuralLddt(structures, 15f),
                  new AlignedPairCoverage(proteins)));
      AlgorithmFactory factory = factory(three, 50, 1.0);

      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> factory.create(AlgorithmType.MOEA_D));
    }
  }

  @Nested
  @DisplayName("When stopping by time")
  class StoppingByTime {

    @ParameterizedTest
    @EnumSource(
        value = AlgorithmType.class,
        names = {"NSGA_II", "SMS_EMOA", "RANDOM_SEEDED"})
    @DisplayName("given a time limit, when run, then it runs past the evaluations until the time")
    void givenTimeLimit_whenRun_thenItRunsUntilTheTime(AlgorithmType type) {
      // Arrange
      SearchRun run =
          new AlgorithmFactory(
                  problem,
                  seeds,
                  new AlgorithmConfiguration(
                      EVALUATIONS, POPULATION_SIZE, POPULATION_SIZE, 0.8, 1.0, 1,
                      MutationType.TM_SAFE_GAP_SHIFT, Duration.ofMillis(300), ""),
                  structures,
                  new SeededRandom(5))
              .create(type);

      // Act
      run.run();

      // Assert
      assertTrue(run.runtimeMillis() >= 300, type + ": " + run.runtimeMillis());
      assertTrue(run.evaluations() > EVALUATIONS, type + ": " + run.evaluations());
    }

    @Test
    @DisplayName("given a non-positive time, when configuring, then an exception is thrown")
    void givenNonPositiveTime_whenConfiguring_thenIllegalArgumentExceptionIsThrown() {
      // Act & Assert
      assertThrows(
          IllegalArgumentException.class,
          () ->
              new AlgorithmConfiguration(
                  EVALUATIONS, POPULATION_SIZE, POPULATION_SIZE, 0.8, 1.0, 1,
                  MutationType.TM_SAFE_GAP_SHIFT, Duration.ZERO, ""));
    }
  }

  @Nested
  @DisplayName("When returning an external archive")
  class ExternalArchive {

    private AlgorithmFactory archiving() {
      return new AlgorithmFactory(
          problem,
          seeds,
          new AlgorithmConfiguration(
              EVALUATIONS * 2, POPULATION_SIZE, POPULATION_SIZE, 0.8, 1.0, 1,
              MutationType.TM_SAFE_GAP_SHIFT, null,
              "--algorithmResult externalArchive --populationSizeWithArchive "
                  + POPULATION_SIZE + " --archiveType crowdingDistanceArchive"),
          structures,
          new SeededRandom(5));
    }

    @Test
    @DisplayName("given NSGA-II, when run, then the front is bounded and mutually non-dominated")
    void givenNsgaTwo_whenRun_thenFrontIsBoundedAndMutuallyNonDominated() {
      // Arrange
      SearchRun run = archiving().create(AlgorithmType.NSGA_II);

      // Act
      run.run();

      // Assert
      List<MsaSolution> front = run.front();
      assertFalse(front.isEmpty());
      assertTrue(front.size() <= POPULATION_SIZE, "archive size " + front.size());
      for (MsaSolution first : front) {
        for (MsaSolution second : front) {
          assertFalse(dominates(first.objectives(), second.objectives()));
        }
      }
    }

    @Test
    @DisplayName("given another algorithm, when building it with the archive, then it fails")
    void givenAnotherAlgorithm_whenBuildingWithArchive_thenIllegalArgumentExceptionIsThrown() {
      // Arrange
      AlgorithmFactory factory = archiving();

      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> factory.create(AlgorithmType.SMS_EMOA));
    }

    private static boolean dominates(double[] first, double[] second) {
      boolean better = false;
      for (int index = 0; index < first.length; index++) {
        if (first[index] > second[index]) {
          return false;
        }
        better |= first[index] < second[index];
      }
      return better;
    }
  }
}
