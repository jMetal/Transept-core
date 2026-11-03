package org.tmm2structalign.algorithm;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Set;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.EnumSource;
import org.tmm2structalign.alignment.GapColumnCompactor;
import org.tmm2structalign.alignment.GapIntervals;
import org.tmm2structalign.alignment.Msa;
import org.tmm2structalign.alignment.MsaSolution;
import org.tmm2structalign.model.ProteinSet;
import org.tmm2structalign.objective.TopologyAwareSumOfPairs;
import org.tmm2structalign.operator.TmSafeGapShiftMutation;
import org.tmm2structalign.problem.StructuralMsaProblem;
import org.tmm2structalign.testing.Alignments;
import org.tmm2structalign.testing.SeededRandom;

@DisplayName("Unit tests for class SeededPopulationCreation")
class SeededPopulationCreationTest {

  private static final int SIZE = 30;

  private ProteinSet proteins;
  private StructuralMsaProblem problem;
  private List<Msa> seeds;

  @BeforeEach
  void setUp() {
    proteins =
        Alignments.proteins(
            "ACDEFGHIKLMN:OOMMMMOOMMMO", "ACDEGHIKLMN:OOMMMOOMMMO", "ACEFGHIKLN:OOMMMOOMMO");
    problem =
        new StructuralMsaProblem(
            "p", proteins, List.of(TopologyAwareSumOfPairs.withDefaultParameters()));
    seeds =
        List.of(
            Msa.fromAlignedRows(
                proteins, List.of("ACDEFGHIKLMN--", "ACDEGHIKLMN---", "ACEFGHIKLN----")),
            Msa.fromAlignedRows(
                proteins, List.of("--ACDEFGHIKLMN", "---ACDEGHIKLMN", "----ACEFGHIKLN")),
            Msa.fromAlignedRows(
                proteins, List.of("ACDEFGHIKLMN", "ACDE-GHIKLMN", "AC-EFGHIKL-N")));
  }

  private SeededPopulationCreation creation(InitializationType type, long seed) {
    return new SeededPopulationCreation(
        problem,
        seeds,
        SIZE,
        new SeededRandom(seed),
        type,
        new TmSafeGapShiftMutation(1.0, Alignments.unstructured(proteins), 60,
            new SeededRandom(seed + 100)));
  }

  @Nested
  @DisplayName("When creating the population with each strategy")
  class Strategies {

    @Test
    @DisplayName("given the crossover strategy, when creating, then it matches the original")
    void givenCrossoverStrategy_whenCreating_thenItMatchesTheOriginal() {
      // Arrange
      SeededPopulationCreation original =
          new SeededPopulationCreation(problem, seeds, SIZE, new SeededRandom(3));

      // Act
      List<MsaSolution> expected = original.create();
      List<MsaSolution> actual = creation(InitializationType.CROSSOVER, 3).create();

      // Assert
      assertEquals(
          expected.stream().map(MsaSolution::variables).toList(),
          actual.stream().map(MsaSolution::variables).toList());
    }

    @ParameterizedTest
    @EnumSource(InitializationType.class)
    @DisplayName("given any strategy, when creating, then the seeds come first and size is exact")
    void givenAnyStrategy_whenCreating_thenSeedsComeFirstAndSizeIsExact(InitializationType type) {
      // Act
      List<MsaSolution> population = creation(type, 5).create();

      // Assert
      assertEquals(SIZE, population.size());
      for (int seed = 0; seed < seeds.size(); seed++) {
        assertEquals(seeds.get(seed).rows(), population.get(seed).variables());
      }
      population.forEach(solution -> assertEquals(3, solution.toMsa().size()));
    }

    @Test
    @DisplayName("given the distinct strategy, when creating, then no two members are equivalent")
    void givenDistinctStrategy_whenCreating_thenNoTwoMembersAreEquivalent() {
      // Act
      List<MsaSolution> population = creation(InitializationType.CROSSOVER_DISTINCT, 7).create();

      // Assert: equivalence ignores gap-only columns, which evaluation removes
      Set<List<GapIntervals>> distinct = new HashSet<>();
      for (MsaSolution solution : population) {
        distinct.add(GapColumnCompactor.compact(List.copyOf(solution.variables())));
      }
      assertEquals(SIZE, distinct.size());
    }

    @Test
    @DisplayName("given the perturbed strategy with the TM-safe mutation, then no helix is broken")
    void givenPerturbedStrategyWithTmSafeMutation_whenCreating_thenNoHelixIsBroken() {
      // Arrange: only the seeds padded at the ends, which have no gaps inside helices
      SeededPopulationCreation creation =
          new SeededPopulationCreation(
              problem,
              seeds.subList(0, 2),
              SIZE,
              new SeededRandom(9),
              InitializationType.PERTURBED,
              new TmSafeGapShiftMutation(
                  1.0, Alignments.unstructured(proteins), 60, new SeededRandom(19)));

      // Act
      List<MsaSolution> population = creation.create();

      // Assert
      int breaking = 0;
      for (MsaSolution solution : population) {
        for (int sequence = 0; sequence < proteins.size(); sequence++) {
          breaking +=
              Alignments.helixBreakingRuns(
                  proteins.get(sequence), solution.variables().get(sequence));
        }
      }
      assertEquals(0, breaking);
    }

    @Test
    @DisplayName("given a strategy that mutates and no mutation, when creating, then it fails")
    void givenStrategyThatMutatesAndNoMutation_whenCreating_thenIllegalArgumentExceptionIsThrown() {
      // Act & Assert
      assertThrows(
          IllegalArgumentException.class,
          () ->
              new SeededPopulationCreation(
                  problem, seeds, SIZE, new SeededRandom(1), InitializationType.PERTURBED, null));
    }
  }

  @Nested
  @DisplayName("When spreading the seeds for MOEA/D")
  class Spreading {

    @Test
    @DisplayName("given three seeds in ten members, when spreading, then they sit at 0, 5 and 9")
    void givenThreeSeedsInTenMembers_whenSpreading_thenTheySitAtEvenlySpacedPositions() {
      // Arrange
      List<MsaSolution> population = new ArrayList<>();
      for (int index = 0; index < 10; index++) {
        population.add(Alignments.solution(proteins, seeds.get(index % 3).alignedRow(0),
            seeds.get(index % 3).alignedRow(1), seeds.get(index % 3).alignedRow(2)));
      }

      // Act
      List<MsaSolution> spread = SeededPopulationCreation.spreadSeeds(population, 3);

      // Assert
      assertSame(population.get(0), spread.get(0));
      assertSame(population.get(1), spread.get(5));
      assertSame(population.get(2), spread.get(9));
      assertSame(population.get(3), spread.get(1));
      assertSame(population.get(9), spread.get(8));
      Set<MsaSolution> identities = Collections.newSetFromMap(new IdentityHashMap<>());
      identities.addAll(spread);
      assertEquals(10, identities.size());
    }

    @Test
    @DisplayName("given a single seed, when spreading, then it sits in the middle")
    void givenSingleSeed_whenSpreading_thenItSitsInTheMiddle() {
      // Arrange
      List<MsaSolution> population = new ArrayList<>();
      for (int index = 0; index < 5; index++) {
        population.add(problem.createSolution(seeds.get(index % 3)));
      }

      // Act
      List<MsaSolution> spread = SeededPopulationCreation.spreadSeeds(population, 1);

      // Assert
      assertSame(population.get(0), spread.get(2));
    }
  }
}
