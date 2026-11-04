package org.tmm2structalign.objective;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.tmm2structalign.alignment.Msa;
import org.tmm2structalign.model.ProteinSet;
import org.tmm2structalign.testing.Alignments;

@DisplayName("Unit tests for class AlignedPairCoverage")
class AlignedPairCoverageTest {

  private static final double TOLERANCE = 1e-12;

  private ProteinSet proteins;
  private AlignedPairCoverage coverage;

  @BeforeEach
  void setUp() {
    proteins = Alignments.proteins("ACDEF:OOMMM", "ACDE:OOMM", "ACD:OOM");
    coverage = new AlignedPairCoverage(proteins);
  }

  @Nested
  @DisplayName("When evaluating alignments")
  class Evaluating {

    @Test
    @DisplayName("given every residue of the shorter proteins paired, then coverage is 1")
    void givenEveryShorterResiduePaired_whenEvaluating_thenCoverageIsOne() {
      // Arrange: maxima 4 + 3 + 3 = 10, all reached
      Msa msa = Msa.fromAlignedRows(proteins, List.of("ACDEF", "ACDE-", "ACD--"));

      // Act & Assert
      assertEquals(1.0, coverage.evaluate(msa), TOLERANCE);
    }

    @Test
    @DisplayName("given proteins in disjoint columns, then coverage is 0")
    void givenProteinsInDisjointColumns_whenEvaluating_thenCoverageIsZero() {
      // Arrange
      Msa msa =
          Msa.fromAlignedRows(
              proteins, List.of("ACDEF-------", "-----ACDE---", "---------ACD"));

      // Act & Assert
      assertEquals(0.0, coverage.evaluate(msa), TOLERANCE);
    }

    @Test
    @DisplayName("given a partial overlap, then aligned pairs are divided by the maxima")
    void givenPartialOverlap_whenEvaluating_thenAlignedPairsAreDividedByMaxima() {
      // Arrange: column occupancies 3, 2, 2, 2, 2, 1 -> 3 + 1 + 1 + 1 + 1 = 7 pairs of 10
      Msa msa = Msa.fromAlignedRows(proteins, List.of("ACDEF-", "AC-D-E", "A-C-D-"));

      // Act & Assert
      assertEquals(0.7, coverage.evaluate(msa), TOLERANCE);
    }

    @Test
    @DisplayName("given gap-only columns, then the value does not change")
    void givenGapOnlyColumns_whenEvaluating_thenValueDoesNotChange() {
      // Arrange
      Msa compact = Msa.fromAlignedRows(proteins, List.of("ACDEF", "ACDE-", "ACD--"));
      Msa padded = Msa.fromAlignedRows(proteins, List.of("AC-DEF-", "AC-DE--", "AC-D---"));

      // Act & Assert
      assertEquals(coverage.evaluate(compact), coverage.evaluate(padded), TOLERANCE);
      assertTrue(coverage.isMaximized());
      assertEquals(AlignedPairCoverage.NAME, coverage.name());
    }
  }

  @Nested
  @DisplayName("When creating the objective or using it wrongly")
  class Invalid {

    @Test
    @DisplayName("given a single protein, when creating, then an exception is thrown")
    void givenSingleProtein_whenCreating_thenIllegalArgumentExceptionIsThrown() {
      // Act & Assert
      assertThrows(
          IllegalArgumentException.class,
          () -> new AlignedPairCoverage(Alignments.proteins("ACD:OOM")));
    }

    @Test
    @DisplayName("given an alignment of other proteins, when evaluating, then it fails")
    void givenAlignmentOfOtherProteins_whenEvaluating_thenIllegalArgumentExceptionIsThrown() {
      // Arrange
      ProteinSet others = Alignments.proteins("ACDEF:OOMMM", "ACDE:OOMM", "ACD:OOM");
      Msa msa = Msa.fromAlignedRows(others, List.of("ACDEF", "ACDE-", "ACD--"));

      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> coverage.evaluate(msa));
    }
  }
}
