package org.transept.operator;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

import java.util.List;
import java.util.Random;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.transept.alignment.MsaSolution;
import org.transept.model.ProteinSet;
import org.transept.testing.Alignments;
import org.transept.testing.SeededRandom;

@DisplayName("Unit tests for class MergeAdjacentGapRunsMutation")
class MergeAdjacentGapRunsMutationTest {

  private RandomSource random;

  @BeforeEach
  void setUp() {
    random = mock(RandomSource.class);
  }

  @Nested
  @DisplayName("When merging two runs")
  class Merging {

    @Test
    @DisplayName("given the left run is kept, when merging, then the right run moves next to it")
    void givenLeftRunIsKept_whenMerging_thenRightRunMovesNextToIt() {
      // Arrange: the row has runs at 1 and 4-5; both boundaries are valid, draw 0.2 keeps left
      when(random.nextDouble()).thenReturn(0.0, 0.2);
      when(random.nextInt(0, 0)).thenReturn(0);
      ProteinSet proteins = Alignments.proteins("ACDEF:OOOOO", "ACDEFGHI:OOOOOOOO");
      MsaSolution solution = Alignments.solution(proteins, "A-CD--EF", "ACDEFGHI");

      // Act
      new MergeAdjacentGapRunsMutation(1.0, random).execute(solution);

      // Assert
      assertEquals(List.of("A---CDEF", "ACDEFGHI"), Alignments.rows(solution));
    }

    @Test
    @DisplayName("given the right run is kept, when merging, then the left run moves next to it")
    void givenRightRunIsKept_whenMerging_thenLeftRunMovesNextToIt() {
      // Arrange: draw 0.7 keeps the right boundary
      when(random.nextDouble()).thenReturn(0.0, 0.7);
      when(random.nextInt(0, 0)).thenReturn(0);
      ProteinSet proteins = Alignments.proteins("ACDEF:OOOOO", "ACDEFGHI:OOOOOOOO");
      MsaSolution solution = Alignments.solution(proteins, "A-CD--EF", "ACDEFGHI");

      // Act
      new MergeAdjacentGapRunsMutation(1.0, random).execute(solution);

      // Assert
      assertEquals(List.of("ACD---EF", "ACDEFGHI"), Alignments.rows(solution));
    }

    @Test
    @DisplayName("given a run that breaks a helix, when merging, then that run is moved")
    void givenRunThatBreaksHelix_whenMerging_thenThatRunIsMoved() {
      // Arrange: A-C breaks the helix A C (M M); the run after D is valid, so it stays
      when(random.nextDouble()).thenReturn(0.0);
      when(random.nextInt(0, 0)).thenReturn(0);
      ProteinSet proteins = Alignments.proteins("ACDEF:MMOOO", "ACDEFGHI:OOOOOOOO");
      MsaSolution solution = Alignments.solution(proteins, "A-CD--EF", "ACDEFGHI");

      // Act
      new MergeAdjacentGapRunsMutation(1.0, random).execute(solution);

      // Assert
      assertEquals(List.of("ACD---EF", "ACDEFGHI"), Alignments.rows(solution));
      assertEquals(0, Alignments.helixBreakingRuns(proteins.get(0), solution.variables().get(0)));
    }
  }

  @Test
  @DisplayName("given random alignments, when merging, then no helix is broken")
  void givenRandomAlignments_whenMerging_thenNoHelixIsBroken() {
    // Arrange
    ProteinSet proteins =
        Alignments.proteins("ACDEFGHIKLMNPQRS:OOMMMMMMOOMMMMMO", "ACDEFGHIKLMN:OMMMMMOOMMMO");
    Random generator = new Random(11);
    MergeAdjacentGapRunsMutation mutation =
        new MergeAdjacentGapRunsMutation(1.0, new SeededRandom(5));

    for (int iteration = 0; iteration < 500; iteration++) {
      MsaSolution solution = Alignments.randomSolution(proteins, 10, generator);
      final int before = breaking(proteins, solution);
      final int runs = runs(solution);
      final int length = solution.alignmentLength();

      // Act
      mutation.execute(solution);

      // Assert
      assertTrue(breaking(proteins, solution) <= before, "iteration " + iteration);
      assertEquals(length, solution.alignmentLength());
      assertTrue(runs(solution) == runs - 1 || runs(solution) == runs, "iteration " + iteration);
      assertEquals(2, solution.toMsa().size());
    }
  }

  @Test
  @DisplayName("given no row with two runs, when mutating, then it returns without changes")
  void givenNoRowWithTwoRuns_whenMutating_thenItReturnsWithoutChanges() {
    // Arrange: the original operator looped forever here
    when(random.nextDouble()).thenReturn(0.0);
    ProteinSet proteins = Alignments.proteins("ACD:OOO", "ACD:OOO");
    MsaSolution solution = Alignments.solution(proteins, "AC-D", "ACD-");

    // Act
    new MergeAdjacentGapRunsMutation(1.0, random).execute(solution);

    // Assert
    assertEquals(List.of("AC-D", "ACD-"), Alignments.rows(solution));
  }

  @Test
  @DisplayName("given invalid arguments, when creating or mutating, then it fails")
  void givenInvalidArguments_whenCreatingOrMutating_thenIllegalArgumentExceptionIsThrown() {
    assertThrows(IllegalArgumentException.class, () -> new MergeAdjacentGapRunsMutation(-0.1));
    MergeAdjacentGapRunsMutation mutation = new MergeAdjacentGapRunsMutation(0.3);
    assertThrows(IllegalArgumentException.class, () -> mutation.execute(null));
    assertEquals(0.3, mutation.mutationProbability());
  }

  private static int breaking(ProteinSet proteins, MsaSolution solution) {
    int breaking = 0;
    for (int sequence = 0; sequence < proteins.size(); sequence++) {
      breaking +=
          Alignments.helixBreakingRuns(proteins.get(sequence), solution.variables().get(sequence));
    }
    return breaking;
  }

  private static int runs(MsaSolution solution) {
    return solution.variables().stream().mapToInt(row -> row.runCount()).sum();
  }
}
