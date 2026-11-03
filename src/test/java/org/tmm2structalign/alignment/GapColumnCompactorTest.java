package org.tmm2structalign.alignment;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.util.List;
import java.util.stream.Stream;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

@DisplayName("Unit tests for class GapColumnCompactor")
class GapColumnCompactorTest {

  @Nested
  @DisplayName("When removing gap-only columns")
  class Compaction {

    static Stream<Arguments> alignments() {
      return Stream.of(
          // Review finding A2, two sequences: the original code removed only 2 of 4 columns.
          Arguments.of(
              List.of("----ACDEFGHIKL", "----ACDEFGHIKL"),
              List.of("ACDEFGHIKL", "ACDEFGHIKL")),
          // Review finding A2, three sequences with several runs in the last row.
          Arguments.of(
              List.of("------ACDEFGHIKL", "------ACDEFGHIKL", "--A--CDEFGHIKL--"),
              List.of("--ACDEFGHIKL", "--ACDEFGHIKL", "ACDEFGHIKL--")),
          Arguments.of(List.of("A--C", "A--C"), List.of("AC", "AC")),
          Arguments.of(List.of("AC--", "A-C-"), List.of("AC-", "A-C")),
          Arguments.of(List.of("A-C", "-BC"), List.of("A-C", "-BC")));
    }

    @ParameterizedTest
    @MethodSource("alignments")
    @DisplayName("given an alignment, when compacting, then only gap-only columns are removed")
    void givenAlignment_whenCompacting_thenOnlyGapOnlyColumnsAreRemoved(
        List<String> aligned, List<String> expected) {
      // Arrange
      List<GapIntervals> rows = encode(aligned);

      // Act
      List<GapIntervals> compacted = GapColumnCompactor.compact(rows);

      // Assert
      assertEquals(encode(expected), compacted);
      assertEquals(0, GapColumnCompactor.countGapOnlyColumns(compacted));
    }

    @Test
    @DisplayName("given an alignment without gap-only columns, when compacting, then it is kept")
    void givenNoGapOnlyColumns_whenCompacting_thenInputIsReturned() {
      // Arrange
      List<GapIntervals> rows = encode(List.of("A-C", "AB-"));

      // Act
      List<GapIntervals> compacted = GapColumnCompactor.compact(rows);

      // Assert
      assertSame(rows, compacted);
    }

    @Test
    @DisplayName("given an alignment, when counting, then gap-only columns are counted")
    void givenAlignment_whenCounting_thenGapOnlyColumnsAreCounted() {
      // Act & Assert
      assertEquals(3, GapColumnCompactor.countGapOnlyColumns(encode(List.of("-A--", "-B--"))));
    }
  }

  @Nested
  @DisplayName("When the input is invalid")
  class InvalidInput {

    @Test
    @DisplayName("given no rows, when compacting, then an exception is thrown")
    void givenNoRows_whenCompacting_thenIllegalArgumentExceptionIsThrown() {
      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> GapColumnCompactor.compact(List.of()));
    }

    @Test
    @DisplayName("given rows of different lengths, when compacting, then an exception is thrown")
    void givenDifferentLengths_whenCompacting_thenIllegalArgumentExceptionIsThrown() {
      // Arrange
      List<GapIntervals> rows = encode(List.of("A-", "A--"));

      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> GapColumnCompactor.compact(rows));
    }
  }

  private static List<GapIntervals> encode(List<String> aligned) {
    return aligned.stream().map(GapIntervals::fromAlignedRow).toList();
  }
}
