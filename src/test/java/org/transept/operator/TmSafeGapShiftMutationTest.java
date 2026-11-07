package org.transept.operator;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertSame;
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
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.transept.alignment.MsaSolution;
import org.transept.model.ProteinSet;
import org.transept.testing.Alignments;
import org.transept.testing.SeededRandom;

@DisplayName("Unit tests for class TmSafeGapShiftMutation")
class TmSafeGapShiftMutationTest {

  private RandomSource random;

  @BeforeEach
  void setUp() {
    random = mock(RandomSource.class);
  }

  @Nested
  @DisplayName("When protecting transmembrane helices")
  class HelixProtection {

    @ParameterizedTest
    @ValueSource(longs = {1, 2, 3, 4, 5})
    @DisplayName("given random alignments, when mutating, then no run is moved into a helix")
    void givenRandomAlignments_whenMutating_thenNoRunIsMovedIntoHelix(long seed) {
      // Arrange: random alignments usually contain gap runs inside helices
      Random generator = new Random(seed);
      ProteinSet proteins =
          Alignments.proteins(
              "ACDEFGHIKLMNPQRS:OOMMMMMMOOMMMMMO", "ACDEFGHIKLMN:OMMMMMOOMMMO");
      TmSafeGapShiftMutation mutation =
          new TmSafeGapShiftMutation(
              1.0, Alignments.unstructured(proteins), 60, new SeededRandom(seed));

      for (int iteration = 0; iteration < 500; iteration++) {
        MsaSolution solution = Alignments.randomSolution(proteins, 10, generator);
        int before = 0;
        for (int sequence = 0; sequence < 2; sequence++) {
          before +=
              Alignments.helixBreakingRuns(
                  proteins.get(sequence), solution.variables().get(sequence));
        }
        final int length = solution.alignmentLength();

        // Act
        mutation.execute(solution);

        // Assert: moving a run can repair a helix but never break one
        int after = 0;
        for (int sequence = 0; sequence < 2; sequence++) {
          after +=
              Alignments.helixBreakingRuns(
                  proteins.get(sequence), solution.variables().get(sequence));
        }
        assertTrue(after <= before, "iteration " + iteration + ": " + before + " -> " + after);
        assertEquals(length, solution.alignmentLength());
        assertEquals(2, solution.toMsa().size());
      }
    }

    @Test
    @DisplayName("given a protein without structure, when proposing a helix, then it is never used")
    void givenProteinWithoutStructure_whenProposingHelix_thenItIsNeverUsed() {
      // Arrange: valid positions of "ACDEF" (O M M M O) after removing the run at 0: 1, 4 and 5.
      // The mock always returns the last valid position, so the run must end before F.
      when(random.nextDouble()).thenReturn(0.0);
      when(random.nextInt(0, 0)).thenReturn(0);
      when(random.nextInt(0, 2)).thenReturn(2);
      ProteinSet proteins = Alignments.proteins("ACDEF:OMMMO", "ACDEFGH:OOOOOOO");
      TmSafeGapShiftMutation mutation =
          new TmSafeGapShiftMutation(1.0, Alignments.unstructured(proteins), 60, random);
      MsaSolution solution = Alignments.solution(proteins, "--ACDEF", "ACDEFGH");

      // Act
      mutation.execute(solution);

      // Assert
      assertEquals(List.of("ACDEF--", "ACDEFGH"), Alignments.rows(solution));
    }

    @Test
    @DisplayName("given the start of a helix, when proposing it, then the boundary is valid")
    void givenStartOfHelix_whenProposingIt_thenBoundaryIsValid() {
      // Arrange: valid positions are 1 (between A and C, O|M), 4 (E|F) and 5 (end); choose 1
      when(random.nextDouble()).thenReturn(0.0);
      when(random.nextInt(0, 0)).thenReturn(0);
      when(random.nextInt(0, 2)).thenReturn(0);
      ProteinSet proteins = Alignments.proteins("ACDEF:OMMMO", "ACDEFGH:OOOOOOO");
      TmSafeGapShiftMutation mutation =
          new TmSafeGapShiftMutation(1.0, Alignments.fullyStructured(proteins), 60, random);
      MsaSolution solution = Alignments.solution(proteins, "--ACDEF", "ACDEFGH");

      // Act
      mutation.execute(solution);

      // Assert
      assertEquals(List.of("A--CDEF", "ACDEFGH"), Alignments.rows(solution));
    }
  }

  @Nested
  @DisplayName("When choosing what to mutate")
  class Choice {

    @Test
    @DisplayName("given only single gaps, when mutating, then a single gap is moved")
    void givenOnlySingleGaps_whenMutating_thenSingleGapIsMoved() {
      // Arrange: the original operator required a run of two or more gaps
      when(random.nextDouble()).thenReturn(0.0);
      when(random.nextInt(0, 0)).thenReturn(0);
      when(random.nextInt(0, 4)).thenReturn(3);
      ProteinSet proteins = Alignments.proteins("ACDEF:OOOOO", "ACDEFG:OOOOOO");
      TmSafeGapShiftMutation mutation =
          new TmSafeGapShiftMutation(1.0, Alignments.fullyStructured(proteins), 60, random);
      MsaSolution solution = Alignments.solution(proteins, "A-CDEF", "ACDEFG");

      // Act
      mutation.execute(solution);

      // Assert: valid positions without the original (1) are 0, 2, 3, 4 and 5; index 3 is 4
      assertEquals(List.of("ACDE-F", "ACDEFG"), Alignments.rows(solution));
    }

    @Test
    @DisplayName("given no gaps, when mutating, then nothing changes")
    void givenNoGaps_whenMutating_thenNothingChanges() {
      // Arrange
      when(random.nextDouble()).thenReturn(0.0);
      ProteinSet proteins = Alignments.proteins("ACD:OOO", "ACD:OOO");
      TmSafeGapShiftMutation mutation =
          new TmSafeGapShiftMutation(1.0, Alignments.fullyStructured(proteins), 60, random);
      MsaSolution solution = Alignments.solution(proteins, "ACD", "ACD");

      // Act
      MsaSolution mutated = mutation.execute(solution);

      // Assert
      assertSame(solution, mutated);
      assertEquals(List.of("ACD", "ACD"), Alignments.rows(solution));
    }

    @Test
    @DisplayName("given only rejections, when mutating, then a valid position is chosen uniformly")
    void givenOnlyRejections_whenMutating_thenValidPositionIsChosenUniformly() {
      // Arrange: acceptance draws of 0.99 reject every attempt; the fallback picks index 1
      when(random.nextDouble()).thenReturn(0.0, 0.99);
      when(random.nextInt(0, 0)).thenReturn(0);
      when(random.nextInt(0, 2)).thenReturn(0, 1);
      ProteinSet proteins = Alignments.proteins("ACDEF:OMMMO", "ACDEFGH:OOOOOOO");
      TmSafeGapShiftMutation mutation =
          new TmSafeGapShiftMutation(1.0, Alignments.fullyStructured(proteins), 1, random);
      MsaSolution solution = Alignments.solution(proteins, "--ACDEF", "ACDEFGH");

      // Act
      mutation.execute(solution);

      // Assert: valid positions 1, 4, 5; the fallback chose 4 (between E and F)
      assertEquals(List.of("ACDE--F", "ACDEFGH"), Alignments.rows(solution));
    }

    @Test
    @DisplayName("given a failed probability draw, when mutating, then nothing changes")
    void givenFailedProbabilityDraw_whenMutating_thenNothingChanges() {
      // Arrange
      when(random.nextDouble()).thenReturn(0.9);
      ProteinSet proteins = Alignments.proteins("ACD:OOO", "AC:OO");
      TmSafeGapShiftMutation mutation =
          new TmSafeGapShiftMutation(0.5, Alignments.fullyStructured(proteins), 60, random);
      MsaSolution solution = Alignments.solution(proteins, "ACD", "A-C");

      // Act
      mutation.execute(solution);

      // Assert
      assertEquals(List.of("ACD", "A-C"), Alignments.rows(solution));
      assertEquals(0.5, mutation.mutationProbability());
    }
  }

  @Nested
  @DisplayName("When limiting how far a run can move")
  class Reach {

    private static final String RESIDUES = "A".repeat(30);
    private static final String ALIGNED = "A".repeat(10) + "--" + "A".repeat(20);

    private ProteinSet proteins() {
      return Alignments.proteins(
          RESIDUES + ":" + "O".repeat(30), "A".repeat(32) + ":" + "O".repeat(32));
    }

    @ParameterizedTest
    @ValueSource(longs = {1, 2, 3, 4, 5, 6, 7, 8})
    @DisplayName("given a limit, when mutating, then the run moves within it")
    void givenLimit_whenMutating_thenRunMovesWithinIt(long seed) {
      // Arrange
      ProteinSet proteins = proteins();
      TmSafeGapShiftMutation mutation =
          new TmSafeGapShiftMutation(
              1.0, Alignments.fullyStructured(proteins), 60, 3, new SeededRandom(seed));
      MsaSolution solution = Alignments.solution(proteins, ALIGNED, "A".repeat(32));

      // Act
      mutation.execute(solution);

      // Assert: it moved, and at most 3 columns
      int start = Alignments.rows(solution).get(0).indexOf('-');
      assertTrue(start != 10 && Math.abs(start - 10) <= 3, "moved to " + start);
    }

    @Test
    @DisplayName("given no valid position within reach, when mutating, then nothing changes")
    void givenNoValidPositionWithinReach_whenMutating_thenNothingChanges() {
      // Arrange: the only other valid position is at the far end of the helix
      ProteinSet proteins =
          Alignments.proteins(
              "A".repeat(10) + ":" + "M".repeat(10), "A".repeat(11) + ":" + "O".repeat(11));
      TmSafeGapShiftMutation mutation =
          new TmSafeGapShiftMutation(
              1.0, Alignments.fullyStructured(proteins), 60, 3, new SeededRandom(1));
      MsaSolution solution = Alignments.solution(proteins, "-" + "A".repeat(10), "A".repeat(11));

      // Act
      mutation.execute(solution);

      // Assert
      assertEquals("-" + "A".repeat(10), Alignments.rows(solution).get(0));
    }

    @ParameterizedTest
    @ValueSource(longs = {1, 2, 3})
    @DisplayName("given no limit, when mutating, then it behaves as without the parameter")
    void givenNoLimit_whenMutating_thenItBehavesAsWithoutTheParameter(long seed) {
      // Arrange
      ProteinSet proteins = proteins();
      var structures = Alignments.fullyStructured(proteins);
      MsaSolution first = Alignments.solution(proteins, ALIGNED, "A".repeat(32));
      MsaSolution second = Alignments.solution(proteins, ALIGNED, "A".repeat(32));

      // Act
      new TmSafeGapShiftMutation(1.0, structures, 60, new SeededRandom(seed)).execute(first);
      new TmSafeGapShiftMutation(
              1.0,
              structures,
              60,
              TmSafeGapShiftMutation.UNLIMITED_SHIFT,
              new SeededRandom(seed))
          .execute(second);

      // Assert
      assertEquals(Alignments.rows(first), Alignments.rows(second));
    }
  }

  @Nested
  @DisplayName("When the operator is misused")
  class Misuse {

    @Test
    @DisplayName("given invalid arguments, when creating or executing, then exceptions are thrown")
    void givenInvalidArguments_whenCreatingOrExecuting_thenExceptionsAreThrown() {
      // Arrange
      ProteinSet proteins = Alignments.proteins("ACD:OOO", "AC:OO");
      TmSafeGapShiftMutation mutation =
          new TmSafeGapShiftMutation(1.0, Alignments.fullyStructured(proteins));

      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> mutation.execute(null));
      assertThrows(
          IllegalArgumentException.class,
          () -> new TmSafeGapShiftMutation(1.5, Alignments.fullyStructured(proteins)));
      assertThrows(
          IllegalArgumentException.class,
          () -> new TmSafeGapShiftMutation(1.0, Alignments.fullyStructured(proteins), 0, random));
      assertThrows(
          IllegalArgumentException.class,
          () ->
              new TmSafeGapShiftMutation(
                  1.0, Alignments.fullyStructured(proteins), 60, -1, random));
    }
  }
}
