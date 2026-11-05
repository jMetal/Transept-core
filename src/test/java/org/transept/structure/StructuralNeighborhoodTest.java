package org.transept.structure;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

@DisplayName("Unit tests for class StructuralNeighborhood")
class StructuralNeighborhoodTest {

  private static final float[][] DISTANCES = {
    {0f, 4f, 10f, Float.NaN},
    {4f, 0f, 6f, 9f},
    {10f, 6f, 0f, 3f},
    {Float.NaN, 9f, 3f, 0f}
  };

  @Nested
  @DisplayName("When building neighbourhoods")
  class Building {

    @Test
    @DisplayName("given a cutoff, when building, then only strictly closer neighbours are kept")
    void givenCutoff_whenBuilding_thenOnlyStrictlyCloserNeighboursAreKept() {
      // Act
      StructuralNeighborhood neighborhood = StructuralNeighborhood.build(DISTANCES, 9.0f);

      // Assert
      assertEquals(4, neighborhood.residueCount());
      assertArrayEquals(new int[] {1}, neighborhood.neighborIndices(0));
      assertArrayEquals(new int[] {0, 2}, neighborhood.neighborIndices(1));
      assertArrayEquals(new float[] {4f, 6f}, neighborhood.neighborDistances(1));
      assertArrayEquals(new float[] {1f - 4f / 9f, 1f - 6f / 9f}, neighborhood.neighborWeights(1));
      assertArrayEquals(new int[] {2}, neighborhood.neighborIndices(3));
    }
  }

  @Nested
  @DisplayName("When the input is invalid")
  class InvalidInput {

    @ParameterizedTest
    @ValueSource(floats = {0f, -1f, Float.NaN, Float.POSITIVE_INFINITY})
    @DisplayName("given an invalid cutoff, when building, then an exception is thrown")
    void givenInvalidCutoff_whenBuilding_thenIllegalArgumentExceptionIsThrown(float cutoff) {
      // Act & Assert
      assertThrows(
          IllegalArgumentException.class, () -> StructuralNeighborhood.build(DISTANCES, cutoff));
    }

    @Test
    @DisplayName("given a non-square matrix, when building, then an exception is thrown")
    void givenNonSquareMatrix_whenBuilding_thenIllegalArgumentExceptionIsThrown() {
      // Arrange
      float[][] distances = {{0f, 1f}, {1f}};

      // Act & Assert
      assertThrows(
          IllegalArgumentException.class, () -> StructuralNeighborhood.build(distances, 5f));
    }
  }
}
