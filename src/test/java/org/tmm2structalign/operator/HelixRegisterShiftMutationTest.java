package org.tmm2structalign.operator;

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
import org.tmm2structalign.alignment.GapIntervals;
import org.tmm2structalign.alignment.MsaSolution;
import org.tmm2structalign.model.ProteinSequence;
import org.tmm2structalign.model.ProteinSet;
import org.tmm2structalign.testing.Alignments;
import org.tmm2structalign.testing.SeededRandom;

@DisplayName("Unit tests for class HelixRegisterShiftMutation")
class HelixRegisterShiftMutationTest {

  /** Helix DEFGHI (residues 2..7) and a protein without helices. */
  private static final ProteinSet PROTEINS =
      Alignments.proteins("ACDEFGHIKL:OOMMMMMMOO", "ACDEFGHIKL:OOOOOOOOOO");

  private RandomSource random;

  @BeforeEach
  void setUp() {
    random = mock(RandomSource.class);
  }

  @Nested
  @DisplayName("When moving a helix")
  class Moving {

    @Test
    @DisplayName("given gaps only on the left, when mutating, then the helix moves left")
    void givenGapsOnlyOnTheLeft_whenMutating_thenHelixMovesLeft() {
      // Arrange: the only candidate; the shift is 2, the number of gaps on the left
      when(random.nextDouble()).thenReturn(0.0);
      when(random.nextInt(0, 0)).thenReturn(0);
      when(random.nextInt(1, 2)).thenReturn(2);
      MsaSolution solution = Alignments.solution(PROTEINS, "AC--DEFGHIKL", "ACDEFGHIKL--");

      // Act
      new HelixRegisterShiftMutation(1.0, 3, random).execute(solution);

      // Assert: the two gaps move to the right of the helix
      assertEquals(List.of("ACDEFGHI--KL", "ACDEFGHIKL--"), Alignments.rows(solution));
    }

    @Test
    @DisplayName("given gaps on both sides, when moving right, then the right gap moves left")
    void givenGapsOnBothSides_whenMovingRight_thenRightGapMovesToTheLeft() {
      // Arrange: draw 0.9 chooses the right; one gap there, so the shift is 1
      when(random.nextDouble()).thenReturn(0.0, 0.9);
      when(random.nextInt(0, 0)).thenReturn(0);
      when(random.nextInt(1, 1)).thenReturn(1);
      MsaSolution solution = Alignments.solution(PROTEINS, "AC-DEFGHI-KL", "ACDEFGHIKL--");

      // Act
      new HelixRegisterShiftMutation(1.0, 3, random).execute(solution);

      // Assert
      assertEquals(List.of("AC--DEFGHIKL", "ACDEFGHIKL--"), Alignments.rows(solution));
    }

    @Test
    @DisplayName("given a helix with a gap inside, when mutating, then it is not moved")
    void givenHelixWithGapInside_whenMutating_thenItIsNotMoved() {
      // Arrange
      when(random.nextDouble()).thenReturn(0.0);
      MsaSolution solution = Alignments.solution(PROTEINS, "ACDE-FGHIKL-", "ACDEFGHIKL--");

      // Act
      new HelixRegisterShiftMutation(1.0, 3, random).execute(solution);

      // Assert
      assertEquals(List.of("ACDE-FGHIKL-", "ACDEFGHIKL--"), Alignments.rows(solution));
    }
  }

  @Test
  @DisplayName("given two helices, when finding candidates, then each loop ends at the other")
  void givenTwoHelices_whenFindingCandidates_thenEachLoopEndsAtTheOtherHelix() {
    // Arrange: helices ACD and HIK; columns -ACD-EFG--HIK-
    ProteinSequence protein = new ProteinSequence("p", "ACDEFGHIK", "MMMOOOMMM");
    GapIntervals row = GapIntervals.fromAlignedRow("-ACD-EFG--HIK-");

    // Act
    List<HelixRegisterShiftMutation.Candidate> candidates =
        HelixRegisterShiftMutation.candidates(0, protein, row);

    // Assert
    assertEquals(2, candidates.size());
    assertEquals(List.of(0), candidates.get(0).leftGaps());
    assertEquals(List.of(4, 8, 9), candidates.get(0).rightGaps());
    assertEquals(List.of(9, 8, 4), candidates.get(1).leftGaps());
    assertEquals(List.of(13), candidates.get(1).rightGaps());
  }

  @Test
  @DisplayName("given random alignments, when mutating, then helices stay whole and rows valid")
  void givenRandomAlignments_whenMutating_thenHelicesStayWholeAndRowsValid() {
    // Arrange
    ProteinSet proteins =
        Alignments.proteins(
            "ACDEFGHIKLMNPQRSTVWY:OOMMMMMOOOOMMMMMMOOO", "ACDEFGHIKLMNPQ:OMMMMOOOOMMMMO",
            "ACDEFGHIKLMNPQRS:OOOOOOOOOOOOOOOO");
    Random generator = new Random(4);
    HelixRegisterShiftMutation mutation =
        new HelixRegisterShiftMutation(1.0, 3, new SeededRandom(6));

    for (int iteration = 0; iteration < 500; iteration++) {
      MsaSolution solution = Alignments.randomSolution(proteins, 8, generator);
      List<GapIntervals> before = List.copyOf(solution.variables());
      final int breaking = breaking(proteins, solution);
      final int length = solution.alignmentLength();

      // Act
      mutation.execute(solution);

      // Assert: at most one row changes, no helix breaks, the length is kept
      int changed = 0;
      for (int sequence = 0; sequence < proteins.size(); sequence++) {
        changed += before.get(sequence).equals(solution.variables().get(sequence)) ? 0 : 1;
      }
      assertTrue(changed <= 1, "iteration " + iteration);
      assertTrue(breaking(proteins, solution) <= breaking, "iteration " + iteration);
      assertEquals(length, solution.alignmentLength());
      assertEquals(3, solution.toMsa().size());
    }
  }

  @Test
  @DisplayName("given invalid arguments, when creating or mutating, then it fails")
  void givenInvalidArguments_whenCreatingOrMutating_thenIllegalArgumentExceptionIsThrown() {
    assertThrows(
        IllegalArgumentException.class, () -> new HelixRegisterShiftMutation(1.5, 3, random));
    assertThrows(
        IllegalArgumentException.class, () -> new HelixRegisterShiftMutation(1.0, 0, random));
    HelixRegisterShiftMutation mutation = new HelixRegisterShiftMutation(0.2);
    assertThrows(IllegalArgumentException.class, () -> mutation.execute(null));
    assertEquals(0.2, mutation.mutationProbability());
    assertEquals(HelixRegisterShiftMutation.DEFAULT_MAX_SHIFT, mutation.maxShift());
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
