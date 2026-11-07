package org.transept.algorithm;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.ValueSource;
import org.transept.alignment.Msa;
import org.transept.alignment.MsaSolution;
import org.transept.model.ProteinSet;
import org.transept.objective.TopologyAwareSumOfPairs;
import org.transept.operator.TmSafeGapShiftMutation;
import org.transept.problem.StructuralMsaProblem;
import org.transept.testing.Alignments;
import org.transept.testing.SeededRandom;

@DisplayName("Unit tests for SeededPopulationCreation and AlgorithmConfiguration")
class AlgorithmTest {

  private ProteinSet proteins;
  private StructuralMsaProblem problem;
  private List<Msa> seeds;

  @BeforeEach
  void setUp() {
    proteins =
        Alignments.proteins("ACDEFGHIK:OOMMMMOOO", "ACDEGHIK:OOMMMOOO", "ACEFGHIK:OOMMMOOO");
    problem =
        new StructuralMsaProblem(
            "test", proteins, List.of(TopologyAwareSumOfPairs.withDefaultParameters()));
    seeds =
        List.of(
            Msa.fromAlignedRows(proteins, List.of("ACDEFGHIK-", "ACDE-GHIK-", "AC-EFGHIK-")),
            Msa.fromAlignedRows(proteins, List.of("ACDEFGHIK", "ACD-EGHIK", "A-CEFGHIK")),
            Msa.fromAlignedRows(proteins, List.of("-ACDEFGHIK", "-ACDEGHIK-", "--ACEFGHIK")));
  }

  private TmSafeGapShiftMutation mutation() {
    return new TmSafeGapShiftMutation(
        1.0, Alignments.unstructured(proteins), 60, new SeededRandom(99));
  }

  @Nested
  @DisplayName("When creating the initial population")
  class InitialPopulation {

    @ParameterizedTest
    @ValueSource(ints = {1, 2, 3, 4, 7, 10})
    @DisplayName("given a population size, when creating, then exactly that many are created")
    void givenPopulationSize_whenCreating_thenExactlyThatManyAreCreated(int size) {
      // Arrange
      SeededPopulationCreation creation =
          new SeededPopulationCreation(problem, seeds, size, new SeededRandom(size), mutation());

      // Act
      List<MsaSolution> population = creation.create();

      // Assert
      assertEquals(size, population.size());
    }

    @Test
    @DisplayName("given duplicated seeds, when creating, then each distinct seed is used once")
    void givenDuplicatedSeeds_whenCreating_thenEachDistinctSeedIsUsedOnce() {
      // Arrange
      List<Msa> duplicated = List.of(seeds.get(0), seeds.get(1), seeds.get(0));
      SeededPopulationCreation creation =
          new SeededPopulationCreation(problem, duplicated, 2, new SeededRandom(1), mutation());

      // Act
      List<MsaSolution> population = creation.create();

      // Assert
      assertEquals(2, creation.distinctSeeds().size());
      assertEquals(seeds.get(0).rows(), population.get(0).variables());
      assertEquals(seeds.get(1).rows(), population.get(1).variables());
    }

    @Test
    @DisplayName("given a single distinct seed and a larger population, when creating, then fail")
    void givenSingleSeedAndLargerPopulation_whenCreating_thenIllegalArgumentExceptionIsThrown() {
      // Arrange
      SeededPopulationCreation creation =
          new SeededPopulationCreation(
              problem, List.of(seeds.get(0), seeds.get(0)), 4, new SeededRandom(1), mutation());

      // Act & Assert
      assertThrows(IllegalArgumentException.class, creation::create);
    }

    @Test
    @DisplayName("given invalid arguments, when creating the builder, then they are rejected")
    void givenInvalidArguments_whenCreatingBuilder_thenIllegalArgumentExceptionIsThrown() {
      // Act & Assert
      assertThrows(
          IllegalArgumentException.class,
          () -> new SeededPopulationCreation(problem, seeds, 0, new SeededRandom(1), mutation()));
      assertThrows(
          IllegalArgumentException.class,
          () ->
              new SeededPopulationCreation(
                  problem, List.of(), 3, new SeededRandom(1), mutation()));
      assertThrows(
          IllegalArgumentException.class,
          () -> new SeededPopulationCreation(problem, seeds, 3, new SeededRandom(1), null));
    }
  }

  @Nested
  @DisplayName("When configuring the algorithm")
  class Configuration {

    @Test
    @DisplayName("given defaults, when creating, then the published probabilities are used")
    void givenDefaults_whenCreating_thenPublishedProbabilitiesAreUsed() {
      // Act
      AlgorithmConfiguration configuration = AlgorithmConfiguration.withDefaults(1000, 100, 4);

      // Assert
      assertEquals(0.8, configuration.crossoverProbability());
      assertEquals(0.2, configuration.mutationProbability());
      assertEquals(100, configuration.offspringPopulationSize());
    }

    @ParameterizedTest
    @CsvSource({
      "0, 10, 10, 0.8, 0.2, 1",
      "10, 0, 10, 0.8, 0.2, 1",
      "10, 10, 0, 0.8, 0.2, 1",
      "10, 10, 10, 1.5, 0.2, 1",
      "10, 10, 10, 0.8, -0.1, 1",
      "10, 10, 10, 0.8, 0.2, 0"
    })
    @DisplayName("given an invalid parameter, when creating, then an exception is thrown")
    void givenInvalidParameter_whenCreating_thenIllegalArgumentExceptionIsThrown(
        int evaluations, int population, int offspring, double pc, double pm, int threads) {
      // Act & Assert
      assertThrows(
          IllegalArgumentException.class,
          () -> new AlgorithmConfiguration(evaluations, population, offspring, pc, pm, threads));
    }

    @ParameterizedTest
    @CsvSource({"-1, 1", "0, 0"})
    @DisplayName("given invalid operator settings, when creating, then an exception is thrown")
    void givenInvalidOperatorSettings_whenCreating_thenIllegalArgumentExceptionIsThrown(
        int maxShift, int cutCandidates) {
      // Act & Assert
      assertThrows(
          IllegalArgumentException.class,
          () ->
              new AlgorithmConfiguration(
                  10,
                  10,
                  10,
                  0.8,
                  0.2,
                  1,
                  MutationType.TM_SAFE_GAP_SHIFT,
                  null,
                  "",
                  maxShift,
                  cutCandidates));
    }
  }
}
