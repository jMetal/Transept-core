package org.tmm2structalign.operator;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.tmm2structalign.alignment.MsaSolution;
import org.tmm2structalign.model.ProteinSet;
import org.tmm2structalign.testing.Alignments;
import org.tmm2structalign.testing.SeededRandom;

@DisplayName("Unit tests for class TmSafeMultiPointCrossover")
class TmSafeMultiPointCrossoverTest {

  private static final ProteinSet PLAIN = Alignments.proteins("ACDEF:OOOOO", "ACDEF:OOOOO");

  @Nested
  @DisplayName("When recombining")
  class Recombining {

    @Test
    @DisplayName("given two cuts, when crossing, then the middle segment comes from the second")
    void givenTwoCuts_whenCrossing_thenMiddleSegmentComesFromSecondParent() {
      // Arrange: first cut at column 2 (index 1 of 1..5), second at column 5 (index 2 of 3..7)
      RandomSource random = mock(RandomSource.class);
      when(random.nextDouble()).thenReturn(0.0);
      when(random.nextInt(0, 4)).thenReturn(1, 2);
      MsaSolution first = Alignments.solution(PLAIN, "ACDEF--", "A--CDEF");
      MsaSolution second = Alignments.solution(PLAIN, "--ACDEF", "ACDEF--");

      // Act
      List<MsaSolution> children =
          new TmSafeMultiPointCrossover(1.0, 2, random).execute(List.of(first, second));

      // Assert: ACD | ---- from the second | EF from the first, and A-- | CDE | F
      assertEquals(List.of("ACD---EF--", "A--CDE---F"), Alignments.rows(children.get(0)));
    }

    @Test
    @DisplayName("given explicit cuts, when recombining, then it matches the random path")
    void givenExplicitCuts_whenRecombining_thenItMatchesTheRandomPath() {
      // Arrange
      MsaSolution first = Alignments.solution(PLAIN, "ACDEF--", "A--CDEF");
      MsaSolution second = Alignments.solution(PLAIN, "--ACDEF", "ACDEF--");

      // Act
      MsaSolution child =
          new TmSafeMultiPointCrossover(1.0, 2).recombine(first, second, List.of(2, 5));

      // Assert
      assertEquals(List.of("ACD---EF--", "A--CDE---F"), Alignments.rows(child));
    }

    @ParameterizedTest
    @ValueSource(longs = {1, 2, 3, 4})
    @DisplayName("given one point, when crossing, then it equals the single-point crossover")
    void givenOnePoint_whenCrossing_thenItEqualsTheSinglePointCrossover(long seed) {
      // Arrange
      Random generator = new Random(seed);
      MsaSolution first = Alignments.randomSolution(PLAIN, 4, generator);
      MsaSolution second = Alignments.randomSolution(PLAIN, 4, generator);

      // Act
      List<MsaSolution> multi =
          new TmSafeMultiPointCrossover(1.0, 1, new SeededRandom(seed))
              .execute(List.of(first, second));
      List<MsaSolution> single =
          new TmSafeSinglePointCrossover(1.0, new SeededRandom(seed))
              .execute(List.of(first, second));

      // Assert
      for (int child = 0; child < 2; child++) {
        assertEquals(single.get(child).variables(), multi.get(child).variables());
      }
    }

    @Test
    @DisplayName("given no valid cut, when crossing, then the child copies the first parent")
    void givenNoValidCut_whenCrossing_thenChildCopiesTheFirstParent() {
      // Arrange: every inner column of the first parent splits the helix
      ProteinSet helix = Alignments.proteins("ACDE:MMMM", "ACDE:MMMM");
      MsaSolution first = Alignments.solution(helix, "-ACDE", "-ACDE");
      MsaSolution second = Alignments.solution(helix, "ACDE-", "ACDE-");

      // Act
      List<MsaSolution> children =
          new TmSafeMultiPointCrossover(1.0, 3, new SeededRandom(1))
              .execute(List.of(first, second));

      // Assert
      assertEquals(first.variables(), children.get(0).variables());
    }
  }

  @Test
  @DisplayName("given parents without broken helices, when crossing, then children have none")
  void givenParentsWithoutBrokenHelices_whenCrossing_thenChildrenHaveNone() {
    // Arrange
    ProteinSet proteins =
        Alignments.proteins("ACDEFGHIKLMNPQRS:OOMMMMMMOOMMMMMO", "ACDEFGHIKLMN:OMMMMMOOMMMO");
    Random generator = new Random(3);
    List<MsaSolution> parents = new ArrayList<>();
    while (parents.size() < 40) {
      MsaSolution candidate = Alignments.randomSolution(proteins, 8, generator);
      if (breaking(proteins, candidate) == 0) {
        parents.add(candidate);
      }
    }
    TmSafeMultiPointCrossover crossover =
        new TmSafeMultiPointCrossover(1.0, 3, new SeededRandom(9));

    for (int pair = 0; pair + 1 < parents.size(); pair += 2) {
      // Act
      List<MsaSolution> children =
          crossover.execute(List.of(parents.get(pair), parents.get(pair + 1)));

      // Assert
      for (MsaSolution child : children) {
        assertEquals(0, breaking(proteins, child), "pair " + pair);
        assertEquals(2, child.toMsa().size());
      }
    }
  }

  @Test
  @DisplayName("given invalid arguments, when creating or crossing, then it fails")
  void givenInvalidArguments_whenCreatingOrCrossing_thenIllegalArgumentExceptionIsThrown() {
    assertThrows(IllegalArgumentException.class, () -> new TmSafeMultiPointCrossover(1.0, 0));
    assertThrows(IllegalArgumentException.class, () -> new TmSafeMultiPointCrossover(1.1, 2));
    TmSafeMultiPointCrossover crossover = new TmSafeMultiPointCrossover(0.9, 2);
    MsaSolution parent = Alignments.solution(PLAIN, "ACDEF", "ACDEF");
    MsaSolution other =
        Alignments.solution(Alignments.proteins("ACDEF:OOOOO", "ACDEF:OOOOO"), "ACDEF", "ACDEF");
    assertThrows(IllegalArgumentException.class, () -> crossover.execute(List.of(parent)));
    assertThrows(IllegalArgumentException.class, () -> crossover.execute(List.of(parent, other)));
    assertEquals(0.9, crossover.crossoverProbability());
    assertEquals(2, crossover.points());
    assertEquals(2, crossover.numberOfRequiredParents());
    assertEquals(2, crossover.numberOfGeneratedChildren());
  }

  private static int breaking(ProteinSet proteins, MsaSolution solution) {
    int breaking = 0;
    for (int sequence = 0; sequence < proteins.size(); sequence++) {
      breaking +=
          Alignments.helixBreakingRuns(proteins.get(sequence), solution.variables().get(sequence));
    }
    return breaking;
  }
}
