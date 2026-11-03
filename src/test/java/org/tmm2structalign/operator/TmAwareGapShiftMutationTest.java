package org.tmm2structalign.operator;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.ArgumentMatchers.anyInt;
import static org.mockito.ArgumentMatchers.eq;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.never;
import static org.mockito.Mockito.verify;
import static org.mockito.Mockito.when;

import java.util.List;
import java.util.Random;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.tmm2structalign.alignment.MsaSolution;
import org.tmm2structalign.model.ProteinSet;
import org.tmm2structalign.testing.Alignments;
import org.tmm2structalign.testing.SeededRandom;

@DisplayName("Unit tests for class TmAwareGapShiftMutation")
class TmAwareGapShiftMutationTest {

  private ProteinSet proteins;
  private RandomSource random;

  @BeforeEach
  void setUp() {
    proteins = Alignments.proteins("ACDEF:OMMMO", "ACDEFGH:OOOOOOO");
    random = mock(RandomSource.class);
  }

  @Nested
  @DisplayName("When moving a gap block with scripted random numbers")
  class Scripted {

    @Test
    @DisplayName("given a TM proposal then a loop proposal, when mutating, then the loop wins")
    void givenTmThenLoopProposal_whenMutating_thenBlockMovesToTheLoop() {
      // Arrange: apply; sequence 0; run 0; propose D (TM, reject), then F (loop, accept).
      when(random.nextDouble()).thenReturn(0.0, 0.5, 0.5);
      when(random.nextInt(0, 0)).thenReturn(0, 0);
      when(random.nextInt(0, 5)).thenReturn(2, 4);
      TmAwareGapShiftMutation mutation =
          new TmAwareGapShiftMutation(1.0, Alignments.fullyStructured(proteins), 60, random);
      MsaSolution solution = Alignments.solution(proteins, "--ACDEF", "ACDEFGH");

      // Act
      MsaSolution mutated = mutation.execute(solution);

      // Assert
      assertSame(solution, mutated);
      assertEquals(List.of("ACDE--F", "ACDEFGH"), Alignments.rows(mutated));
    }

    @Test
    @DisplayName("given no structure, when proposing a TM residue, then it is usually accepted")
    void givenNoStructure_whenProposingTmResidue_thenItIsAccepted() {
      // Arrange: apply; sequence 0; run 0; propose D with draw 0.94 < 0.95.
      when(random.nextDouble()).thenReturn(0.0, 0.94);
      when(random.nextInt(0, 0)).thenReturn(0, 0);
      when(random.nextInt(0, 5)).thenReturn(2);
      TmAwareGapShiftMutation mutation =
          new TmAwareGapShiftMutation(1.0, Alignments.unstructured(proteins), 60, random);
      MsaSolution solution = Alignments.solution(proteins, "--ACDEF", "ACDEFGH");

      // Act
      mutation.execute(solution);

      // Assert
      assertEquals(List.of("AC--DEF", "ACDEFGH"), Alignments.rows(solution));
    }

    @Test
    @DisplayName("given only rejections, when mutating, then a uniform position is used")
    void givenOnlyRejections_whenMutating_thenUniformFallbackIsUsed() {
      // Arrange: one try, proposal equals the forbidden start (0), then fallback 0, then 5.
      when(random.nextDouble()).thenReturn(0.0);
      when(random.nextInt(0, 0)).thenReturn(0, 0);
      when(random.nextInt(0, 5)).thenReturn(0, 0, 5);
      TmAwareGapShiftMutation mutation =
          new TmAwareGapShiftMutation(1.0, Alignments.fullyStructured(proteins), 1, random);
      MsaSolution solution = Alignments.solution(proteins, "--ACDEF", "ACDEFGH");

      // Act
      mutation.execute(solution);

      // Assert
      assertEquals(List.of("ACDEF--", "ACDEFGH"), Alignments.rows(solution));
    }

    @Test
    @DisplayName("given only single gaps, when mutating, then nothing changes and no loop occurs")
    void givenOnlySingleGaps_whenMutating_thenNothingChanges() {
      // Arrange: review finding M1, the original operator looped forever here.
      when(random.nextDouble()).thenReturn(0.0);
      TmAwareGapShiftMutation mutation =
          new TmAwareGapShiftMutation(1.0, Alignments.fullyStructured(proteins), 60, random);
      MsaSolution solution = Alignments.solution(proteins, "A-CDE-F", "ACDEFGH");

      // Act
      mutation.execute(solution);

      // Assert
      assertEquals(List.of("A-CDE-F", "ACDEFGH"), Alignments.rows(solution));
      verify(random, never()).nextInt(anyInt(), eq(5));
    }

    @Test
    @DisplayName("given a failed probability draw, when mutating, then nothing changes")
    void givenFailedProbabilityDraw_whenMutating_thenNothingChanges() {
      // Arrange
      when(random.nextDouble()).thenReturn(0.5);
      TmAwareGapShiftMutation mutation =
          new TmAwareGapShiftMutation(0.2, Alignments.fullyStructured(proteins), 60, random);
      MsaSolution solution = Alignments.solution(proteins, "--ACDEF", "ACDEFGH");

      // Act
      mutation.execute(solution);

      // Assert
      assertEquals(List.of("--ACDEF", "ACDEFGH"), Alignments.rows(solution));
      assertEquals(0.2, mutation.mutationProbability());
    }
  }

  @Nested
  @DisplayName("When mutating random alignments")
  class Invariants {

    @ParameterizedTest
    @ValueSource(longs = {11, 12, 13, 14, 15})
    @DisplayName("given random alignments, when mutating, then length and residues are kept")
    void givenRandomAlignments_whenMutating_thenLengthAndResiduesAreKept(long seed) {
      // Arrange
      Random generator = new Random(seed);
      ProteinSet longer =
          Alignments.proteins("ACDEFGHIKLMN:OOMMMMOOMMMO", "ACDEFGHIK:OMMMMOOMO");
      TmAwareGapShiftMutation mutation =
          new TmAwareGapShiftMutation(
              1.0, Alignments.fullyStructured(longer), 60, new SeededRandom(seed));

      for (int iteration = 0; iteration < 200; iteration++) {
        MsaSolution solution = Alignments.randomSolution(longer, 8, generator);
        int length = solution.alignmentLength();

        // Act
        mutation.execute(solution);

        // Assert: toMsa() re-validates that every row still spells its protein.
        assertEquals(length, solution.alignmentLength());
        assertEquals(2, solution.toMsa().size());
      }
    }

    @Test
    @DisplayName("given many mutations, when counting targets, then TM residues are avoided")
    void givenManyMutations_whenCountingTargets_thenTmResiduesAreAvoided() {
      // Arrange: residues 0-4 are TM, 5-9 are not; one block of two gaps at the start.
      ProteinSet halfTm = Alignments.proteins("ACDEFGHIKL:MMMMMOOOOO", "ACDEFGHIKLMN:OOOOOOOOOOOO");
      TmAwareGapShiftMutation mutation =
          new TmAwareGapShiftMutation(
              1.0, Alignments.fullyStructured(halfTm), 60, new SeededRandom(42));
      int insideHelix = 0;
      int outsideHelix = 0;

      for (int iteration = 0; iteration < 2000; iteration++) {
        MsaSolution solution = Alignments.solution(halfTm, "--ACDEFGHIKL", "ACDEFGHIKLMN");

        // Act
        mutation.execute(solution);

        // Assert (aggregated below): where does the block start now?
        int start = solution.variables().get(0).start(0);
        if (start >= 1 && start <= 4) {
          insideHelix++;
        } else if (start >= 6 && start <= 9) {
          outsideHelix++;
        }
      }
      assertTrue(outsideHelix > 5 * insideHelix, outsideHelix + " vs " + insideHelix);
    }
  }

  @Nested
  @DisplayName("When the operator is misused")
  class Misuse {

    @Test
    @DisplayName("given a null solution, when mutating, then an exception is thrown")
    void givenNullSolution_whenMutating_thenIllegalArgumentExceptionIsThrown() {
      // Arrange
      TmAwareGapShiftMutation mutation =
          new TmAwareGapShiftMutation(1.0, Alignments.fullyStructured(proteins));

      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> mutation.execute(null));
    }

    @ParameterizedTest
    @ValueSource(doubles = {-0.1, 1.5, Double.NaN})
    @DisplayName("given an invalid probability, when creating, then an exception is thrown")
    void givenInvalidProbability_whenCreating_thenIllegalArgumentExceptionIsThrown(double value) {
      // Act & Assert
      assertThrows(
          IllegalArgumentException.class,
          () -> new TmAwareGapShiftMutation(value, Alignments.fullyStructured(proteins)));
    }

    @Test
    @DisplayName("given zero tries, when creating, then an exception is thrown")
    void givenZeroTries_whenCreating_thenIllegalArgumentExceptionIsThrown() {
      // Act & Assert
      assertThrows(
          IllegalArgumentException.class,
          () ->
              new TmAwareGapShiftMutation(
                  1.0, Alignments.fullyStructured(proteins), 0, random));
    }
  }
}
