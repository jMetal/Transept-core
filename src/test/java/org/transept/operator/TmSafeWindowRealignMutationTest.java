package org.transept.operator;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import java.util.Random;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.transept.alignment.GapIntervals;
import org.transept.alignment.MsaSolution;
import org.transept.model.ProteinSet;
import org.transept.testing.Alignments;
import org.transept.testing.SeededRandom;

@DisplayName("Unit tests for class TmSafeWindowRealignMutation")
class TmSafeWindowRealignMutationTest {

  private static List<String> realign(
      ProteinSet proteins, int sequence, int first, int columns, String... rows) {
    MsaSolution solution = Alignments.solution(proteins, rows);
    GapIntervals row = TmSafeWindowRealignMutation.realign(solution, sequence, first, columns);
    solution.variables().set(sequence, row);
    return Alignments.rows(solution);
  }

  @Test
  @DisplayName("given a misplaced gap, when realigning, then the row matches the others")
  void givenMisplacedGap_whenRealigning_thenRowMatchesTheOthers() {
    // Arrange
    ProteinSet proteins = Alignments.proteins("WWCC:OOOO", "WWCC:OOOO", "WWCC:OOOO");

    // Act
    List<String> rows = realign(proteins, 0, 0, 5, "W-WCC", "WWCC-", "WWCC-");

    // Assert: the other rows do not change
    assertEquals(List.of("WWCC-", "WWCC-", "WWCC-"), rows);
  }

  @Test
  @DisplayName("given a window, when realigning, then nothing changes outside it")
  void givenWindow_whenRealigning_thenNothingChangesOutsideIt() {
    // Arrange
    ProteinSet proteins = Alignments.proteins("AWWCC:OOOOO", "AWWCC:OOOOO");

    // Act: the gap at the start stays, and the one inside the window moves
    List<String> rows = realign(proteins, 0, 2, 5, "-AW-WCC", "A-WWCC-");

    // Assert
    assertEquals("-AWWCC-", rows.get(0));
  }

  @Test
  @DisplayName("given a better score across a helix, when realigning, then the helix stays whole")
  void givenBetterScoreAcrossHelix_whenRealigning_thenHelixStaysWhole() {
    // Arrange: the other row would pull a gap between the two transmembrane residues
    ProteinSet proteins = Alignments.proteins("WWCC:OMMO", "WWCC:OOOO");
    MsaSolution solution = Alignments.solution(proteins, "WWCC-", "WW-CC");

    // Act
    GapIntervals row = TmSafeWindowRealignMutation.realign(solution, 0, 0, 5);

    // Assert
    assertEquals(0, Alignments.helixBreakingRuns(proteins.get(0), row));
    assertEquals(5, row.alignedLength());
  }

  @Test
  @DisplayName("given a window inside a helix, when realigning, then nothing changes")
  void givenWindowInsideHelix_whenRealigning_thenNothingChanges() {
    // Arrange: an inherited gap, with no valid place for it inside the window
    ProteinSet proteins = Alignments.proteins("AAAA:MMMM", "AAAAA:OOOOO");

    // Act
    List<String> rows = realign(proteins, 0, 1, 3, "AA-AA", "AAAAA");

    // Assert
    assertEquals("AA-AA", rows.get(0));
  }

  @ParameterizedTest
  @ValueSource(longs = {1, 2, 3, 4, 5, 6})
  @DisplayName("given random alignments, when mutating, then it stays valid")
  void givenRandomAlignments_whenMutating_thenItIsValidAndBreaksNoMoreHelices(long seed) {
    // Arrange
    ProteinSet proteins =
        Alignments.proteins(
            "ACDEFGHIKLMNPQRS:OOMMMMMMOOMMMMMO",
            "ACDEFGHIKLMN:OMMMMMOOMMMO",
            "ACDEFGHIK:OOOOOOOOO");
    MsaSolution solution = Alignments.randomSolution(proteins, 6, new Random(seed));
    int length = solution.alignmentLength();
    int before = 0;
    for (int row = 0; row < proteins.size(); row++) {
      before += Alignments.helixBreakingRuns(proteins.get(row), solution.variables().get(row));
    }

    // Act
    new TmSafeWindowRealignMutation(1.0, 8, new SeededRandom(seed)).execute(solution);

    // Assert
    int after = 0;
    for (int row = 0; row < proteins.size(); row++) {
      GapIntervals gaps = solution.variables().get(row);
      after += Alignments.helixBreakingRuns(proteins.get(row), gaps);
      assertEquals(proteins.get(row).length(), gaps.residueCount());
      assertEquals(length, gaps.alignedLength());
    }
    assertTrue(after <= before, before + " -> " + after);
  }

  @Test
  @DisplayName("given invalid arguments, when creating or executing, then exceptions are thrown")
  void givenInvalidArguments_whenCreatingOrExecuting_thenExceptionsAreThrown() {
    // Arrange
    SeededRandom random = new SeededRandom(1);

    // Act & Assert
    assertThrows(
        IllegalArgumentException.class, () -> new TmSafeWindowRealignMutation(1.5, 10, random));
    assertThrows(
        IllegalArgumentException.class, () -> new TmSafeWindowRealignMutation(0.5, 1, random));
    assertThrows(
        IllegalArgumentException.class,
        () -> new TmSafeWindowRealignMutation(0.5, 10, random).execute(null));
  }
}
