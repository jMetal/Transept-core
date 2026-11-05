package org.transept.algorithm;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertSame;

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
import org.transept.alignment.GapColumnCompactor;
import org.transept.alignment.GapIntervals;
import org.transept.alignment.Msa;
import org.transept.alignment.MsaSolution;
import org.transept.model.ProteinSet;
import org.transept.objective.TopologyAwareSumOfPairs;
import org.transept.operator.TmSafeGapShiftMutation;
import org.transept.problem.StructuralMsaProblem;
import org.transept.testing.Alignments;
import org.transept.testing.SeededRandom;

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

  private SeededPopulationCreation creation(List<Msa> chosen, long seed) {
    return new SeededPopulationCreation(
        problem,
        chosen,
        SIZE,
        new SeededRandom(seed),
        new TmSafeGapShiftMutation(1.0, Alignments.unstructured(proteins), 60,
            new SeededRandom(seed + 100)));
  }

  @Nested
  @DisplayName("When creating the population")
  class Creating {

    @Test
    @DisplayName("given three seeds, when creating, then the seeds come first and size is exact")
    void givenThreeSeeds_whenCreating_thenSeedsComeFirstAndSizeIsExact() {
      // Act
      List<MsaSolution> population = creation(seeds, 5).create();

      // Assert
      assertEquals(SIZE, population.size());
      for (int seed = 0; seed < seeds.size(); seed++) {
        assertEquals(seeds.get(seed).rows(), population.get(seed).variables());
      }
      population.forEach(solution -> assertEquals(3, solution.toMsa().size()));
    }

    @Test
    @DisplayName("given three seeds, when creating, then no two members are equivalent")
    void givenThreeSeeds_whenCreating_thenNoTwoMembersAreEquivalent() {
      // Act
      List<MsaSolution> population = creation(seeds, 7).create();

      // Assert: equivalence ignores gap-only columns, which evaluation removes
      Set<List<GapIntervals>> distinct = new HashSet<>();
      for (MsaSolution solution : population) {
        distinct.add(GapColumnCompactor.compact(List.copyOf(solution.variables())));
      }
      assertEquals(SIZE, distinct.size());
    }

    @Test
    @DisplayName("given seeds without gaps in helices, when creating, then no helix is broken")
    void givenSeedsWithoutGapsInHelices_whenCreating_thenNoHelixIsBroken() {
      // Arrange: only the seeds padded at the ends, which have no gaps inside helices
      SeededPopulationCreation creation = creation(seeds.subList(0, 2), 9);

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
