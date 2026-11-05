package org.transept.alignment;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.transept.model.ProteinSequence;
import org.transept.model.ProteinSet;

@DisplayName("Unit tests for class MsaSolution")
class MsaSolutionTest {

  private ProteinSet proteins;
  private MsaSolution solution;

  @BeforeEach
  void setUp() {
    proteins =
        new ProteinSet(
            List.of(
                ProteinSequence.withoutTopology("p1", "ACD"),
                ProteinSequence.withoutTopology("p2", "AD")));
    solution =
        new MsaSolution(
            proteins,
            List.of(GapIntervals.fromAlignedRow("-ACD-"), GapIntervals.fromAlignedRow("-A-D-")),
            2);
  }

  @Nested
  @DisplayName("When creating and copying solutions")
  class CreationAndCopy {

    @Test
    @DisplayName("given valid rows, when creating, then the alignment is decoded")
    void givenValidRows_whenCreating_thenAlignmentIsDecoded() {
      // Act
      Msa msa = solution.toMsa();

      // Assert
      assertEquals(5, solution.alignmentLength());
      assertEquals("-A-D-", msa.alignedRow(1));
      assertSame(proteins, solution.proteins());
    }

    @Test
    @DisplayName("given rows that do not fit the proteins, when creating, then it is rejected")
    void givenInconsistentRows_whenCreating_thenIllegalArgumentExceptionIsThrown() {
      // Arrange
      List<GapIntervals> rows = List.of(GapIntervals.none(3), GapIntervals.none(3));

      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> new MsaSolution(proteins, rows, 2));
    }

    @Test
    @DisplayName("given an evaluated solution, when copying, then rows and objectives are copied")
    void givenEvaluatedSolution_whenCopying_thenRowsAndObjectivesAreCopied() {
      // Arrange
      solution.objectives()[0] = 1.5;
      solution.objectives()[1] = -2.0;
      solution.attributes().put("key", "value");

      // Act
      MsaSolution copy = solution.copy();
      copy.objectives()[0] = 0.0;

      // Assert
      assertNotSame(solution, copy);
      assertEquals(solution.variables(), copy.variables());
      assertArrayEquals(new double[] {1.5, -2.0}, solution.objectives());
      assertEquals(-2.0, copy.objectives()[1]);
      assertEquals("value", copy.attributes().get("key"));
    }
  }

  @Nested
  @DisplayName("When compacting gap columns")
  class Compaction {

    @Test
    @DisplayName("given gap-only columns, when compacting, then they are removed from all rows")
    void givenGapOnlyColumns_whenCompacting_thenTheyAreRemoved() {
      // Act
      solution.compactGapColumns();

      // Assert
      assertEquals(3, solution.alignmentLength());
      assertEquals("ACD", solution.toMsa().alignedRow(0));
      assertEquals("A-D", solution.toMsa().alignedRow(1));
    }
  }
}
