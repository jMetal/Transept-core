package org.tmm2structalign.alignment;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.tmm2structalign.model.ProteinSequence;
import org.tmm2structalign.model.ProteinSet;

@DisplayName("Unit tests for class Msa")
class MsaTest {

  private ProteinSet proteins;

  @BeforeEach
  void setUp() {
    proteins =
        new ProteinSet(
            List.of(
                new ProteinSequence("p1", "ACDEF", "OMMMO"),
                new ProteinSequence("p2", "ACDE", "MMMM")));
  }

  @Nested
  @DisplayName("When building an alignment from aligned rows")
  class FromAlignedRows {

    @Test
    @DisplayName("given matching rows, when building, then letters and occupancy are exposed")
    void givenMatchingRows_whenBuilding_thenLettersAndOccupancyAreExposed() {
      // Act
      Msa msa = Msa.fromAlignedRows(proteins, List.of("ACDEF", "AC-DE"));

      // Assert
      assertEquals(2, msa.size());
      assertEquals(5, msa.length());
      assertEquals('-', msa.letter(1, 2));
      assertEquals('D', msa.letter(1, 3));
      assertEquals(3, msa.residueIndex(1, 4));
      assertTrue(msa.isGap(1, 2));
      assertEquals(1, msa.occupancy(2));
      assertEquals("AC-DE", msa.alignedRow(1));
      assertEquals(proteins, msa.proteins());
      assertEquals(2, msa.rows().size());
    }

    @Test
    @DisplayName("given lowercase residues, when building, then case is ignored")
    void givenLowercaseResidues_whenBuilding_thenCaseIsIgnored() {
      // Act
      Msa msa = Msa.fromAlignedRows(proteins, List.of("acdef", "ac-de"));

      // Assert
      assertEquals('A', msa.letter(0, 0));
    }

    @Test
    @DisplayName("given a row with other residues, when building, then an exception is thrown")
    void givenWrongResidues_whenBuilding_thenIllegalArgumentExceptionIsThrown() {
      // Act & Assert
      assertThrows(
          IllegalArgumentException.class,
          () -> Msa.fromAlignedRows(proteins, List.of("ACDEF", "AC-DW")));
    }

    @Test
    @DisplayName("given a missing row, when building, then an exception is thrown")
    void givenMissingRow_whenBuilding_thenIllegalArgumentExceptionIsThrown() {
      // Act & Assert
      assertThrows(
          IllegalArgumentException.class, () -> Msa.fromAlignedRows(proteins, List.of("ACDEF")));
    }
  }

  @Nested
  @DisplayName("When building an alignment from gap layouts")
  class FromGapLayouts {

    @Test
    @DisplayName("given rows of different lengths, when building, then an exception is thrown")
    void givenDifferentLengths_whenBuilding_thenIllegalArgumentExceptionIsThrown() {
      // Arrange
      List<GapIntervals> rows = List.of(GapIntervals.none(5), GapIntervals.of(6, 0, 1));

      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> Msa.of(proteins, rows));
    }

    @Test
    @DisplayName("given a row with the wrong residue count, when building, then it is rejected")
    void givenWrongResidueCount_whenBuilding_thenIllegalArgumentExceptionIsThrown() {
      // Arrange
      List<GapIntervals> rows = List.of(GapIntervals.none(5), GapIntervals.none(5));

      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> Msa.of(proteins, rows));
    }

    @Test
    @DisplayName("given the wrong number of rows, when building, then an exception is thrown")
    void givenWrongNumberOfRows_whenBuilding_thenIllegalArgumentExceptionIsThrown() {
      // Act & Assert
      assertThrows(
          IllegalArgumentException.class, () -> Msa.of(proteins, List.of(GapIntervals.none(5))));
    }
  }

  @Nested
  @DisplayName("When reading transmembrane status")
  class TransmembraneStatus {

    @ParameterizedTest
    @CsvSource({
      // p1 = O M M M O
      "'ACDEF', 0, false",
      "'ACDEF', 1, true",
      // gap between two TM residues (C and D) is TM
      "'AC-DEF', 2, true",
      // gap between O and M residues (A and C) is not TM
      "'A-CDEF', 1, false",
      // leading gap uses the first residue on both sides (O)
      "'-ACDEF', 0, false",
      // trailing gap uses the last residue on both sides (O)
      "'ACDEF-', 5, false"
    })
    @DisplayName("given a column of p1, when reading TM status, then flanking residues decide gaps")
    void givenColumnOfFirstProtein_whenReadingTmStatus_thenFlankingResiduesDecideGaps(
        String alignedFirst, int column, boolean expected) {
      // Arrange
      String alignedSecond = "ACDE" + "-".repeat(alignedFirst.length() - 4);
      Msa msa = Msa.fromAlignedRows(proteins, List.of(alignedFirst, alignedSecond));

      // Act & Assert
      assertEquals(expected, msa.isTransmembrane(0, column));
    }

    @Test
    @DisplayName("given terminal gaps in an all-TM protein, when reading, then the gaps are TM")
    void givenTerminalGapsInTmProtein_whenReading_thenGapsAreTm() {
      // Act
      Msa msa = Msa.fromAlignedRows(proteins, List.of("ACDEF-", "-ACDE-"));

      // Assert
      assertTrue(msa.isTransmembrane(1, 0));
      assertTrue(msa.isTransmembrane(1, 5));
      assertFalse(msa.isTransmembrane(0, 5));
    }
  }
}
