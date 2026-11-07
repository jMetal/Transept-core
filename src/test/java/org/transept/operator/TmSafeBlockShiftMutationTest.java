package org.transept.operator;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.transept.alignment.MsaSolution;
import org.transept.model.ProteinSet;
import org.transept.testing.Alignments;

@DisplayName("Unit tests for class TmSafeBlockShiftMutation")
class TmSafeBlockShiftMutationTest {

  private static final double APPLY_WHOLE_BLOCK = 0.0;
  private static final double LEFT = 0.0;
  private static final double RIGHT = 0.9;

  private RandomSource random;

  @BeforeEach
  void setUp() {
    random = mock(RandomSource.class);
    when(random.nextInt(0, 1)).thenReturn(0);
    when(random.nextInt(0, 0)).thenReturn(0);
    when(random.nextInt(1, 1)).thenReturn(1);
  }

  private MsaSolution mutate(ProteinSet proteins, double direction, String... rows) {
    when(random.nextDouble()).thenReturn(0.0, APPLY_WHOLE_BLOCK, direction);
    MsaSolution solution = Alignments.solution(proteins, rows);
    return new TmSafeBlockShiftMutation(1.0, 1, random).execute(solution);
  }

  @Test
  @DisplayName("given rows sharing a block, when mutating, then they move together")
  void givenRowsSharingBlock_whenMutating_thenTheyMoveTogether() {
    // Arrange
    ProteinSet proteins = Alignments.proteins("ACDE:OOOO", "ACDE:OOOO", "ACDEFG:OOOOOO");

    // Act
    MsaSolution child = mutate(proteins, RIGHT, "AC--DE", "AC--DE", "ACDEFG");

    // Assert: the row with residues in the block stays
    assertEquals(List.of("ACD--E", "ACD--E", "ACDEFG"), Alignments.rows(child));
  }

  @Test
  @DisplayName("given a helix on one side, when mutating, then the block moves the other way")
  void givenHelixOnOneSide_whenMutating_thenBlockMovesTheOtherWay() {
    // Arrange: moving right would separate two transmembrane residues
    ProteinSet proteins = Alignments.proteins("ACDE:OOMM", "ACDEFG:OOOOOO");

    // Act
    MsaSolution child = mutate(proteins, RIGHT, "AC--DE", "ACDEFG");

    // Assert
    assertEquals("A--CDE", Alignments.rows(child).get(0));
  }

  @Test
  @DisplayName("given a row that cannot move, when mutating, then the others still move")
  void givenRowThatCannotMove_whenMutating_thenTheOthersStillMove() {
    // Arrange
    ProteinSet proteins = Alignments.proteins("ACDE:OOOO", "ACDE:OOMM");

    // Act
    MsaSolution child = mutate(proteins, RIGHT, "AC--DE", "AC--DE");

    // Assert
    assertEquals(List.of("ACD--E", "AC--DE"), Alignments.rows(child));
  }

  @Test
  @DisplayName("given no valid direction, when mutating, then nothing changes")
  void givenNoValidDirection_whenMutating_thenNothingChanges() {
    // Arrange
    ProteinSet proteins = Alignments.proteins("ACDE:MMMM", "ACDEF:OOOOO");

    // Act
    MsaSolution child = mutate(proteins, LEFT, "-ACDE", "ACDEF");

    // Assert
    assertEquals(List.of("-ACDE", "ACDEF"), Alignments.rows(child));
  }

  @Test
  @DisplayName("given no gaps or a failed draw, when mutating, then nothing changes")
  void givenNoGapsOrFailedDraw_whenMutating_thenNothingChanges() {
    // Arrange
    ProteinSet withGap = Alignments.proteins("ACDE:OOOO", "ACD:OOO");
    ProteinSet withoutGaps = Alignments.proteins("ACDE:OOOO", "ACDE:OOOO");
    when(random.nextDouble()).thenReturn(1.0, 0.0);
    TmSafeBlockShiftMutation mutation = new TmSafeBlockShiftMutation(0.5, 1, random);

    // Act & Assert: a failed draw, and then a row without gaps
    assertEquals(
        List.of("ACDE", "ACD-"),
        Alignments.rows(mutation.execute(Alignments.solution(withGap, "ACDE", "ACD-"))));
    assertEquals(
        List.of("ACDE", "ACDE"),
        Alignments.rows(mutation.execute(Alignments.solution(withoutGaps, "ACDE", "ACDE"))));
  }

  @Test
  @DisplayName("given invalid arguments, when creating or executing, then exceptions are thrown")
  void givenInvalidArguments_whenCreatingOrExecuting_thenExceptionsAreThrown() {
    // Act & Assert
    assertThrows(
        IllegalArgumentException.class, () -> new TmSafeBlockShiftMutation(1.5, 1, random));
    assertThrows(
        IllegalArgumentException.class, () -> new TmSafeBlockShiftMutation(0.5, 0, random));
    assertThrows(
        IllegalArgumentException.class,
        () -> new TmSafeBlockShiftMutation(0.5, 1, random).execute(null));
  }
}
