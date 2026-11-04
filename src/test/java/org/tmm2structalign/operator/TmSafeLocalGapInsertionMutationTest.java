package org.tmm2structalign.operator;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

import java.util.List;
import java.util.Random;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.tmm2structalign.alignment.MsaSolution;
import org.tmm2structalign.model.ProteinSequence;
import org.tmm2structalign.model.ProteinSet;
import org.tmm2structalign.testing.Alignments;
import org.tmm2structalign.testing.SeededRandom;

@DisplayName("Unit tests for class TmSafeLocalGapInsertionMutation")
class TmSafeLocalGapInsertionMutationTest {

  @Nested
  @DisplayName("When inserting gaps")
  class Inserting {

    @Test
    @DisplayName("given two groups, when inserting, then their gaps are a few columns apart")
    void givenTwoGroups_whenInserting_thenTheirGapsAreFewColumnsApart() {
      // Arrange: position 4, offset +1, the first row alone in the first group
      RandomSource random = mock(RandomSource.class);
      when(random.nextDouble()).thenReturn(0.0, 0.7, 0.2, 0.8, 0.8);
      when(random.nextInt(0, 9)).thenReturn(4);
      when(random.nextInt(1, 3)).thenReturn(1);
      ProteinSet proteins =
          Alignments.proteins("ACDEFGHIK:OOOOOOOOO", "ACDEFGHIK:OOOOOOOOO", "ACDEFGHIK:OOOOOOOOO");
      MsaSolution solution = Alignments.solution(proteins, "ACDEFGHIK", "ACDEFGHIK", "ACDEFGHIK");

      // Act
      new TmSafeLocalGapInsertionMutation(1.0, 1, 3, random).execute(solution);

      // Assert: only F is misaligned
      assertEquals(
          List.of("ACDE-FGHIK", "ACDEF-GHIK", "ACDEF-GHIK"), Alignments.rows(solution));
    }

    @ParameterizedTest
    @CsvSource({"2, 1", "3, 4", "0, 0", "5, 5"})
    @DisplayName("given a target inside a helix, when searching, then the nearest valid is used")
    void givenTargetInsideHelix_whenSearching_thenNearestValidPositionIsUsed(
        int target, int expected) {
      // Arrange: ACDEF (O M M M O); positions 2 and 3 split the helix
      ProteinSequence protein = new ProteinSequence("p", "ACDEF", "OMMMO");

      // Act
      int position = TmSafeLocalGapInsertionMutation.nearestValid(protein, new boolean[5], target);

      // Assert
      assertEquals(expected, position);
    }
  }

  @Test
  @DisplayName("given random alignments, when mutating, then the damage is local and TM-safe")
  void givenRandomAlignments_whenMutating_thenDamageIsLocalAndTmSafe() {
    // Arrange
    ProteinSet plain =
        Alignments.proteins("ACDEFGHIKLMNPQ:OOOOOOOOOOOOOO", "ACDEFGHIKLMN:OOOOOOOOOOOO",
            "ACDEFGHIKL:OOOOOOOOOO");
    ProteinSet helices =
        Alignments.proteins("ACDEFGHIKLMNPQRS:OOMMMMMMOOMMMMMO", "ACDEFGHIKLMN:OMMMMMOOMMMO");
    Random generator = new Random(5);
    int window = 3;
    TmSafeLocalGapInsertionMutation mutation =
        new TmSafeLocalGapInsertionMutation(1.0, 2, window, new SeededRandom(8));

    for (int iteration = 0; iteration < 300; iteration++) {
      // Without helices, two rows that shared a column lose at most `window` of those columns
      MsaSolution solution = Alignments.randomSolution(plain, 6, generator);
      int[][] before = residueColumns(solution);
      mutation.execute(solution);
      int[][] after = residueColumns(solution);
      for (int first = 0; first < 3; first++) {
        for (int second = first + 1; second < 3; second++) {
          assertTrue(
              brokenPairs(before, after, first, second) <= window, "iteration " + iteration);
        }
      }

      // With helices, no new break appears and the alignment grows by gapLength columns
      MsaSolution withHelices = Alignments.randomSolution(helices, 8, generator);
      final int breaking = breaking(helices, withHelices);
      final int length = withHelices.alignmentLength();
      mutation.execute(withHelices);
      assertEquals(breaking, breaking(helices, withHelices), "iteration " + iteration);
      assertEquals(length + 2, withHelices.alignmentLength());
      assertEquals(2, withHelices.toMsa().size());
    }
  }

  @Test
  @DisplayName("given probability zero, when mutating, then nothing changes")
  void givenProbabilityZero_whenMutating_thenNothingChanges() {
    // Arrange
    ProteinSet two = Alignments.proteins("ACD:OOO", "ACD:OOO");
    MsaSolution pair = Alignments.solution(two, "AC-D", "ACD-");

    // Act
    new TmSafeLocalGapInsertionMutation(0.0, 1, 3, new SeededRandom(1)).execute(pair);

    // Assert
    assertEquals(List.of("AC-D", "ACD-"), Alignments.rows(pair));
  }

  @Test
  @DisplayName("given invalid arguments, when creating or mutating, then it fails")
  void givenInvalidArguments_whenCreatingOrMutating_thenIllegalArgumentExceptionIsThrown() {
    RandomSource random = new SeededRandom(1);
    assertThrows(
        IllegalArgumentException.class,
        () -> new TmSafeLocalGapInsertionMutation(1.2, 1, 3, random));
    assertThrows(
        IllegalArgumentException.class,
        () -> new TmSafeLocalGapInsertionMutation(1.0, 0, 3, random));
    assertThrows(
        IllegalArgumentException.class,
        () -> new TmSafeLocalGapInsertionMutation(1.0, 1, 0, random));
    TmSafeLocalGapInsertionMutation mutation = new TmSafeLocalGapInsertionMutation(0.4);
    assertThrows(IllegalArgumentException.class, () -> mutation.execute(null));
    assertEquals(0.4, mutation.mutationProbability());
    assertEquals(1, mutation.gapLength());
    assertEquals(TmSafeLocalGapInsertionMutation.DEFAULT_WINDOW, mutation.window());
  }

  /** Column of every residue of every row. */
  private static int[][] residueColumns(MsaSolution solution) {
    int[][] columns = new int[solution.variables().size()][];
    for (int sequence = 0; sequence < columns.length; sequence++) {
      int[] residueAtColumn = solution.variables().get(sequence).columnToResidue();
      columns[sequence] = new int[solution.proteins().get(sequence).length()];
      for (int column = 0; column < residueAtColumn.length; column++) {
        if (residueAtColumn[column] >= 0) {
          columns[sequence][residueAtColumn[column]] = column;
        }
      }
    }
    return columns;
  }

  /** Residue pairs of two rows aligned in the same column before and not after. */
  private static int brokenPairs(int[][] before, int[][] after, int first, int second) {
    int broken = 0;
    for (int left = 0; left < before[first].length; left++) {
      for (int right = 0; right < before[second].length; right++) {
        if (before[first][left] == before[second][right]
            && after[first][left] != after[second][right]) {
          broken++;
        }
      }
    }
    return broken;
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
