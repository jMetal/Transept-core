package org.tmm2structalign.evolver;

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
import org.tmm2structalign.alignment.Msa;
import org.tmm2structalign.alignment.MsaSolution;
import org.tmm2structalign.model.ProteinSet;
import org.tmm2structalign.objective.StructuralLddt;
import org.tmm2structalign.objective.TopologyAwareSumOfPairs;
import org.tmm2structalign.problem.StructuralMsaProblem;
import org.tmm2structalign.structure.ProteinStructures;
import org.tmm2structalign.testing.Alignments;
import org.tmm2structalign.testing.SeededRandom;
import org.uma.evolver.algorithm.BaseLevelAlgorithm;
import org.uma.jmetal.component.algorithm.EvolutionaryAlgorithm;

@DisplayName("Unit tests for classes MsaSpea2 and MsaGeneticAlgorithm")
class MsaSpea2AndGaTest {

  private static final int POPULATION_SIZE = 8;
  private static final int EVALUATIONS = 48;

  private ProteinSet proteins;
  private ProteinStructures structures;
  private List<Msa> seeds;

  @BeforeEach
  void setUp() {
    proteins =
        Alignments.proteins("ACDEFGHIK:OOMMMMOOO", "ACDEGHIK:OOMMMOOO", "ACEFGHIK:OOMMMOOO");
    structures = Alignments.fullyStructured(proteins);
    seeds =
        List.of(
            Msa.fromAlignedRows(proteins, List.of("ACDEFGHIK-", "ACDE-GHIK-", "AC-EFGHIK-")),
            Msa.fromAlignedRows(proteins, List.of("ACDEFGHIK", "ACD-EGHIK", "A-CEFGHIK")),
            Msa.fromAlignedRows(proteins, List.of("-ACDEFGHIK", "-ACDEGHIK-", "--ACEFGHIK")));
  }

  private StructuralMsaProblem problem(boolean single) {
    StructuralLddt lddt = new StructuralLddt(structures, 15f);
    return new StructuralMsaProblem(
        "test",
        proteins,
        single ? List.of(lddt) : List.of(TopologyAwareSumOfPairs.withDefaultParameters(), lddt));
  }

  private MsaRunContext context(int threads, Duration maxTime) {
    return new MsaRunContext(seeds, structures, new SeededRandom(5), threads, maxTime);
  }

  @SuppressWarnings("unchecked")
  private static EvolutionaryAlgorithm<MsaSolution> build(
      BaseLevelAlgorithm<MsaSolution> template,
      StructuralMsaProblem problem,
      String resource,
      Map<String, String> changes) {
    Map<String, String> values = new LinkedHashMap<>(MsaConfigurations.defaults(resource));
    values.put("offspringPopulationSize", Integer.toString(POPULATION_SIZE));
    values.putAll(changes);
    return (EvolutionaryAlgorithm<MsaSolution>)
        template
            .createInstance(problem, EVALUATIONS)
            .parse(MsaConfigurations.arguments(values))
            .build();
  }

  @Nested
  @DisplayName("When running SPEA2")
  class Spea2 {

    private EvolutionaryAlgorithm<MsaSolution> spea2(Map<String, String> changes, int threads) {
      var template =
          new MsaSpea2(
              POPULATION_SIZE,
              MsaConfigurations.parameterSpace(MsaSpea2.PARAMETER_SPACE),
              context(threads, null));
      return build(template, problem(false), MsaSpea2.PARAMETER_SPACE, changes);
    }

    @ParameterizedTest
    @ValueSource(ints = {1, 2})
    @DisplayName("given the default configuration, when run, then it keeps the population size")
    void givenDefaultConfiguration_whenRun_thenItKeepsThePopulationSize(int threads) {
      // Arrange
      var algorithm = spea2(Map.of(), threads);

      // Act
      algorithm.run();

      // Assert
      assertTrue(algorithm.numberOfEvaluations() >= EVALUATIONS);
      assertEquals(POPULATION_SIZE, algorithm.result().size());
    }

    @Test
    @DisplayName("given NSGA-II's ranking and density, when run, then it completes")
    void givenNsgaTwoRankingAndDensity_whenRun_thenItCompletes() {
      // Arrange
      var algorithm =
          spea2(
              Map.of(
                  "ranking", "dominanceRanking",
                  "densityEstimator", "crowdingDistance",
                  "removalPolicy", "oneShot"),
              1);

      // Act
      algorithm.run();

      // Assert
      assertFalse(algorithm.result().isEmpty());
    }

    @Test
    @DisplayName("given the default configuration, when read, then it is SPEA2")
    void givenDefaultConfiguration_whenRead_thenItIsSpea2() {
      // Act
      Map<String, String> defaults = MsaConfigurations.defaults(MsaSpea2.PARAMETER_SPACE);

      // Assert
      assertEquals("strengthRanking", defaults.get("ranking"));
      assertEquals("knn", defaults.get("densityEstimator"));
      assertEquals("true", defaults.get("knnNormalizeObjectives"));
      assertEquals("sequential", defaults.get("removalPolicy"));
      MsaConfigurations.requireKnownNames(
          defaults, MsaConfigurations.parameterSpace(MsaSpea2.PARAMETER_SPACE));
    }
  }

  @Nested
  @DisplayName("When running the genetic algorithm")
  class GeneticAlgorithm {

    private EvolutionaryAlgorithm<MsaSolution> ga(
        StructuralMsaProblem problem, Map<String, String> changes, Duration maxTime) {
      var template =
          new MsaGeneticAlgorithm(
              POPULATION_SIZE,
              MsaConfigurations.parameterSpace(MsaGeneticAlgorithm.PARAMETER_SPACE),
              context(1, maxTime));
      return build(template, problem, MsaGeneticAlgorithm.PARAMETER_SPACE, changes);
    }

    @Test
    @DisplayName("given (mu, lambda) with more offspring, when run, then it completes")
    void givenMuCommaLambdaWithMoreOffspring_whenRun_thenItCompletes() {
      // Arrange
      var algorithm =
          ga(
              problem(true),
              Map.of("gaReplacement", "muCommaLambda", "offspringPopulationSize", "16"),
              null);

      // Act
      algorithm.run();

      // Assert
      assertEquals(POPULATION_SIZE, algorithm.result().size());
      assertEquals(MsaGeneticAlgorithm.NAME, algorithm.name());
    }

    @Test
    @DisplayName("given a time limit, when run, then it stops by time")
    void givenTimeLimit_whenRun_thenItStopsByTime() {
      // Arrange
      var algorithm = ga(problem(true), Map.of(), Duration.ofMillis(300));

      // Act
      algorithm.run();

      // Assert
      assertTrue(algorithm.totalComputingTime() >= 300);
    }

    @Test
    @DisplayName("given two objectives, when building, then an exception is thrown")
    void givenTwoObjectives_whenBuilding_thenIllegalArgumentExceptionIsThrown() {
      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> ga(problem(false), Map.of(), null));
    }
  }
}
