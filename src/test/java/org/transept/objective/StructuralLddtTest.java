package org.transept.objective;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.transept.alignment.Msa;
import org.transept.model.ProteinSequence;
import org.transept.model.ProteinSet;
import org.transept.structure.MappingStatus;
import org.transept.structure.ProteinStructures;
import org.transept.structure.StructureMapping;

@DisplayName("Unit tests for class StructuralLddt")
class StructuralLddtTest {

  /** Three residues on a line: P(0,1) = 1, P(1,2) = 1, P(0,2) = 2. */
  private static final float[][] STRAIGHT = {{0f, 1f, 2f}, {1f, 0f, 1f}, {2f, 1f, 0f}};

  /** Same, but residues 0 and 2 are 5 Å apart. */
  private static final float[][] BENT = {{0f, 1f, 5f}, {1f, 0f, 1f}, {5f, 1f, 0f}};

  private static final double TOLERANCE = 1e-6;

  @Nested
  @DisplayName("When scoring hand-computed alignments")
  class HandComputed {

    /*
     * STRAIGHT -> BENT. Anchor 0: neighbours 1 (d=1, w=14/15, exact) and 2 (d=2, w=13/15,
     * |2-5|=3 -> 0.25): (14 + 3.25) / 27. Anchor 1: both neighbours exact: 1. Anchor 2: as
     * anchor 0. BENT -> STRAIGHT. Anchor 0: neighbours 1 (w=14/15, exact) and 2 (d=5, w=10/15,
     * 0.25): (14 + 2.5) / 24. Anchor 1: 1. Anchor 2: as anchor 0.
     */
    private static final double STRAIGHT_TO_BENT = 17.25 / 27.0;
    private static final double BENT_TO_STRAIGHT = 16.5 / 24.0;

    @Test
    @DisplayName("given two structures, when scoring, then the mean of the anchors is returned")
    void givenTwoStructures_whenScoring_thenMeanOfAnchorsIsReturned() {
      // Arrange
      StructuralLddt objective = objective(STRAIGHT, BENT);

      // Act
      double value = objective.evaluate(identityAlignment(2, 3));

      // Assert
      double expected = (2 * STRAIGHT_TO_BENT + 1 + 2 * BENT_TO_STRAIGHT + 1) / 6.0;
      assertEquals(expected, value, TOLERANCE);
    }

  }

  @Nested
  @DisplayName("When checking general properties")
  class Properties {

    @Test
    @DisplayName("given identical structures, when aligned identically, then the score is one")
    void givenIdenticalStructures_whenAlignedIdentically_thenScoreIsOne() {
      // Arrange
      StructuralLddt objective = objective(STRAIGHT, STRAIGHT);

      // Act & Assert
      assertEquals(1.0, objective.evaluate(identityAlignment(2, 3)), TOLERANCE);
      assertEquals("lDDT", objective.name());
      assertTrue(objective.isMaximized());
    }

    @Test
    @DisplayName("given an extra gap-only column, when scoring, then the length does not matter")
    void givenExtraGapOnlyColumn_whenScoring_thenLengthDoesNotMatter() {
      // Arrange: the original aggregation divided by the length and gave 3/4 here
      StructuralLddt objective = objective(STRAIGHT, STRAIGHT);
      Msa padded = Msa.fromAlignedRows(proteins(2, 3), List.of("AAA-", "AAA-"));

      // Act & Assert
      assertEquals(1.0, objective.evaluate(padded), TOLERANCE);
    }

    @Test
    @DisplayName("given a protein without structure, when scoring, then it is skipped")
    void givenProteinWithoutStructure_whenScoring_thenItIsSkipped() {
      // Arrange
      StructuralLddt objective = objective(STRAIGHT, STRAIGHT, null);

      // Act
      double value = objective.evaluate(identityAlignment(3, 3));

      // Assert
      assertEquals(1.0, value, TOLERANCE);
    }

    @Test
    @DisplayName("given no usable structures, when scoring, then the score is zero")
    void givenNoUsableStructures_whenScoring_thenScoreIsZero() {
      // Arrange
      StructuralLddt objective = objective(null, null);

      // Act & Assert
      assertEquals(0.0, objective.evaluate(identityAlignment(2, 3)));
    }

    @Test
    @DisplayName("given misaligned residues, when scoring, then unmatched neighbours lower it")
    void givenMisalignedResidues_whenScoring_thenScoreIsLower() {
      // Arrange
      StructuralLddt objective = objective(STRAIGHT, STRAIGHT);
      Msa shifted = Msa.fromAlignedRows(proteins(2, 3), List.of("AAA-", "-AAA"));

      // Act
      double value = objective.evaluate(shifted);

      // Assert
      assertTrue(value < 1.0);
    }

    @Test
    @DisplayName("given a different number of sequences, when scoring, then it is rejected")
    void givenDifferentNumberOfSequences_whenScoring_thenIllegalArgumentExceptionIsThrown() {
      // Arrange
      StructuralLddt objective = objective(STRAIGHT, STRAIGHT);

      // Act & Assert
      assertThrows(
          IllegalArgumentException.class, () -> objective.evaluate(identityAlignment(3, 3)));
    }
  }

  private static StructuralLddt objective(float[][]... matrices) {
    List<StructureMapping> mappings = new ArrayList<>();
    List<float[][]> distances = new ArrayList<>();
    for (int protein = 0; protein < matrices.length; protein++) {
      float[][] matrix = matrices[protein];
      mappings.add(
          matrix == null
              ? new StructureMapping(
                  "p" + protein, MappingStatus.MISSING_FILES, 0, 0, new int[0], "")
              : new StructureMapping(
                  "p" + protein, MappingStatus.NUMBERED, 1, 1, new int[] {0, 1, 2}, ""));
      distances.add(matrix);
    }
    return new StructuralLddt(
        new ProteinStructures(mappings, distances), StructuralLddt.DEFAULT_CUTOFF);
  }

  private static ProteinSet proteins(int count, int length) {
    List<ProteinSequence> proteins = new ArrayList<>();
    for (int protein = 0; protein < count; protein++) {
      proteins.add(ProteinSequence.withoutTopology("p" + protein, "A".repeat(length)));
    }
    return new ProteinSet(proteins);
  }

  private static Msa identityAlignment(int count, int length) {
    List<String> rows = new ArrayList<>();
    for (int protein = 0; protein < count; protein++) {
      rows.add("A".repeat(length));
    }
    return Msa.fromAlignedRows(proteins(count, length), rows);
  }
}
