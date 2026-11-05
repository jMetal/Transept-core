package org.transept.alignment;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.ValueSource;

@DisplayName("Unit tests for class GapIntervals")
class GapIntervalsTest {

  @Nested
  @DisplayName("When encoding and decoding aligned rows")
  class EncodingAndDecoding {

    @ParameterizedTest
    @ValueSource(strings = {"ACDE", "--ACDE", "ACDE--", "A--C-DE", "-A-C-D-E-", "----"})
    @DisplayName("given an aligned row, when encoding and decoding, then the row is restored")
    void givenAlignedRow_whenEncodingAndDecoding_thenRowIsRestored(String alignedRow) {
      // Arrange
      String residues = alignedRow.replace("-", "");

      // Act
      GapIntervals row = GapIntervals.fromAlignedRow(alignedRow);

      // Assert
      assertEquals(alignedRow, row.decode(residues));
      assertEquals(alignedRow.length(), row.alignedLength());
      assertEquals(residues.length(), row.residueCount());
    }

    @Test
    @DisplayName("given consecutive gaps, when encoding, then they form a single run")
    void givenConsecutiveGaps_whenEncoding_thenTheyFormSingleRun() {
      // Act
      GapIntervals row = GapIntervals.fromAlignedRow("A---C-D");

      // Assert
      assertEquals(List.of(1, 3, 5, 5), row.toList());
      assertEquals(2, row.runCount());
      assertEquals(4, row.gapCount());
    }

    @Test
    @DisplayName("given a row, when mapping columns to residues, then gaps map to minus one")
    void givenRow_whenMappingColumnsToResidues_thenGapsMapToMinusOne() {
      // Arrange
      GapIntervals row = GapIntervals.fromAlignedRow("-AC--D");

      // Act
      int[] residueAtColumn = row.columnToResidue();

      // Assert
      assertArrayEquals(new int[] {-1, 0, 1, -1, -1, 2}, residueAtColumn);
    }

    @Test
    @DisplayName("given residues of the wrong length, when decoding, then an exception is thrown")
    void givenWrongResidueCount_whenDecoding_thenIllegalArgumentExceptionIsThrown() {
      // Arrange
      GapIntervals row = GapIntervals.fromAlignedRow("A-C");

      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> row.decode("ACD"));
    }
  }

  @Nested
  @DisplayName("When creating rows from start/end pairs")
  class CreationFromPairs {

    @ParameterizedTest
    @CsvSource({
      "5, '0,5'", // end out of range
      "5, '-1,1'", // negative start
      "5, '3,2'", // start after end
      "6, '0,1,1,3'", // overlapping runs
      "6, '0,1,2,3'", // adjacent runs, must be merged
      "6, '3,4,0,1'", // unsorted runs
      "6, '0'" // incomplete pair
    })
    @DisplayName("given invalid runs, when creating a row, then an exception is thrown")
    void givenInvalidRuns_whenCreatingRow_thenIllegalArgumentExceptionIsThrown(
        int alignedLength, String pairs) {
      // Arrange
      int[] values = parse(pairs);

      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> GapIntervals.of(alignedLength, values));
    }

    @Test
    @DisplayName("given a negative length, when creating a row, then an exception is thrown")
    void givenNegativeLength_whenCreatingRow_thenIllegalArgumentExceptionIsThrown() {
      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> GapIntervals.none(-1));
    }

    @Test
    @DisplayName("given a flat list, when creating a row, then it matches the varargs factory")
    void givenFlatList_whenCreatingRow_thenItMatchesTheVarargsFactory() {
      // Act
      GapIntervals fromList = GapIntervals.fromList(8, List.of(0, 1, 4, 6));

      // Assert
      assertEquals(GapIntervals.of(8, 0, 1, 4, 6), fromList);
      assertEquals(3, fromList.runLength(1));
      assertEquals(4, fromList.start(1));
      assertEquals(6, fromList.end(1));
    }

    @Test
    @DisplayName("given no gaps, when creating a row, then every column is a residue")
    void givenNoGaps_whenCreatingRow_thenEveryColumnIsResidue() {
      // Act
      GapIntervals row = GapIntervals.none(3);

      // Assert
      assertEquals(0, row.runCount());
      assertEquals(3, row.residueCount());
      assertArrayEquals(new boolean[3], row.toMask());
    }
  }

  @Nested
  @DisplayName("When querying gap positions")
  class Queries {

    private final GapIntervals row = GapIntervals.of(10, 0, 1, 4, 4, 8, 9);

    @ParameterizedTest
    @CsvSource({"0,true", "1,true", "2,false", "3,false", "4,true", "5,false", "8,true", "9,true"})
    @DisplayName("given a column, when checking for a gap, then the run membership is returned")
    void givenColumn_whenCheckingForGap_thenRunMembershipIsReturned(int column, boolean gap) {
      // Act & Assert
      assertEquals(gap, row.isGapAt(column));
    }

    @Test
    @DisplayName("given a row, when converting to a mask and back, then the row is unchanged")
    void givenRow_whenConvertingToMaskAndBack_thenRowIsUnchanged() {
      // Act
      GapIntervals roundTrip = GapIntervals.fromMask(row.toMask());

      // Assert
      assertEquals(row, roundTrip);
      assertEquals(row.hashCode(), roundTrip.hashCode());
    }
  }

  @Nested
  @DisplayName("When comparing rows")
  class Equality {

    @Test
    @DisplayName("given rows with different lengths, when comparing, then they are not equal")
    void givenDifferentLengths_whenComparing_thenRowsAreNotEqual() {
      // Act & Assert
      assertNotEquals(GapIntervals.of(5, 0, 0), GapIntervals.of(6, 0, 0));
      assertNotEquals(GapIntervals.of(5, 0, 0), GapIntervals.of(5, 1, 1));
      assertFalse(GapIntervals.none(2).equals((Object) "--"));
    }

    @Test
    @DisplayName("given a row, when printing, then the runs are shown")
    void givenRow_whenPrinting_thenRunsAreShown() {
      // Act & Assert
      assertTrue(GapIntervals.of(5, 1, 2).toString().contains("[1, 2]"));
    }
  }

  private static int[] parse(String pairs) {
    String[] parts = pairs.split(",");
    int[] values = new int[parts.length];
    for (int index = 0; index < parts.length; index++) {
      values[index] = Integer.parseInt(parts[index].trim());
    }
    return values;
  }
}
