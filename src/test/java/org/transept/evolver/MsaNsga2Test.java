package org.transept.evolver;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.time.Duration;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.transept.alignment.Msa;
import org.transept.alignment.MsaSolution;
import org.transept.model.ProteinSet;
import org.transept.objective.StructuralLddt;
import org.transept.objective.TopologyAwareSumOfPairs;
import org.transept.problem.StructuralMsaProblem;
import org.transept.structure.ProteinStructures;
import org.transept.testing.Alignments;
import org.transept.testing.SeededRandom;
import org.uma.jmetal.component.algorithm.EvolutionaryAlgorithm;

@DisplayName("Unit tests for class MsaNsga2")
class MsaNsga2Test {

  private static final int POPULATION_SIZE = 8;
  private static final int EVALUATIONS = 48;

  private ProteinSet proteins;
  private ProteinStructures structures;
  private StructuralMsaProblem problem;
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

  private EvolutionaryAlgorithm<MsaSolution> build(
      Map<String, String> changes, int threads, Duration maxTime) {
    Map<String, String> values =
        new LinkedHashMap<>(MsaConfigurations.defaults(MsaNsga2.PARAMETER_SPACE));
    values.put("offspringPopulationSize", Integer.toString(POPULATION_SIZE));
    values.putAll(changes);
    var template =
        new MsaNsga2(
            POPULATION_SIZE,
            MsaConfigurations.parameterSpace(MsaNsga2.PARAMETER_SPACE),
            new MsaRunContext(seeds, structures, new SeededRandom(5), threads, maxTime));
    return (EvolutionaryAlgorithm<MsaSolution>)
        template
            .createInstance(problem, EVALUATIONS)
            .parse(MsaConfigurations.arguments(values))
            .build();
  }

  @Nested
  @DisplayName("When running configured instances")
  class Running {

    @ParameterizedTest
    @ValueSource(ints = {1, 2})
    @DisplayName("given the default configuration, when run, then it starts from the seeds")
    void givenDefaultConfiguration_whenRun_thenItStartsFromTheSeedsAndUsesTheBudget(int threads) {
      // Arrange
      var algorithm = build(Map.of(), threads, null);

      // Act
      algorithm.run();

      // Assert
      assertTrue(algorithm.numberOfEvaluations() >= EVALUATIONS);
      assertEquals(POPULATION_SIZE, algorithm.result().size());
      algorithm.result().forEach(solution -> assertEquals(3, solution.toMsa().size()));
    }

    @Test
    @DisplayName("given a changed configuration, when run, then it completes")
    void givenChangedConfiguration_whenRun_thenItCompletes() {
      // Arrange
      var algorithm =
          build(
              Map.of(
                  "selection", "random",
                  "mutation", MsaOperators.RANDOM_CHOICE,
                  "mutationProbability", "1.0"),
              1,
              null);

      // Act
      algorithm.run();

      // Assert
      assertFalse(algorithm.result().isEmpty());
    }

    @ParameterizedTest
    @ValueSource(ints = {1, 2})
    @DisplayName("given an external archive, when run, then the archive is bounded")
    void givenExternalArchive_whenRun_thenTheArchiveIsTheBoundedResult(int threads) {
      // Arrange
      var algorithm =
          build(
              Map.of(
                  "algorithmResult", "externalArchive",
                  "populationSizeWithArchive", "10",
                  "archiveType", "crowdingDistanceArchive"),
              threads,
              null);

      // Act
      algorithm.run();

      // Assert: the archive has the population size; the population is populationSizeWithArchive
      assertFalse(algorithm.result().isEmpty());
      assertTrue(algorithm.result().size() <= POPULATION_SIZE);
    }

    @Test
    @DisplayName("given a time limit, when run, then it stops by time instead of evaluations")
    void givenTimeLimit_whenRun_thenItStopsByTime() {
      // Arrange
      var algorithm = build(Map.of(), 1, Duration.ofMillis(300));

      // Act
      algorithm.run();

      // Assert
      assertTrue(algorithm.totalComputingTime() >= 300);
      assertTrue(algorithm.numberOfEvaluations() > EVALUATIONS);
    }
  }

  @Nested
  @DisplayName("When creating the run context")
  class Context {

    @Test
    @DisplayName("given invalid values, when creating the context, then an exception is thrown")
    void givenInvalidValues_whenCreatingContext_thenIllegalArgumentExceptionIsThrown() {
      // Arrange
      var random = new SeededRandom(1);

      // Act & Assert
      assertThrows(
          IllegalArgumentException.class,
          () -> new MsaRunContext(List.of(), structures, random, 1, null));
      assertThrows(
          IllegalArgumentException.class,
          () -> new MsaRunContext(seeds, structures, random, 0, null));
      assertThrows(
          IllegalArgumentException.class,
          () -> new MsaRunContext(seeds, structures, random, 1, Duration.ZERO));
    }
  }
}
