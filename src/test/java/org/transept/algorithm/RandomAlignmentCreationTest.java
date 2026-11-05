package org.transept.algorithm;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.transept.alignment.MsaSolution;
import org.transept.objective.TopologyAwareSumOfPairs;
import org.transept.problem.StructuralMsaProblem;
import org.transept.testing.Alignments;
import org.transept.testing.SeededRandom;

@DisplayName("Unit tests for class RandomAlignmentCreation")
class RandomAlignmentCreationTest {

  private StructuralMsaProblem problem;

  @BeforeEach
  void setUp() {
    problem =
        new StructuralMsaProblem(
            "p",
            Alignments.proteins("ACDEFGHIKL:OOMMMMMMOO", "ACDEFGH:OOMMMOO"),
            List.of(TopologyAwareSumOfPairs.withDefaultParameters()));
  }

  @Nested
  @DisplayName("When creating random alignments")
  class Creation {

    @Test
    @DisplayName("given a batch size, when creating, then that many valid alignments are created")
    void givenBatchSize_whenCreating_thenThatManyValidAlignmentsAreCreated() {
      // Arrange
      RandomAlignmentCreation creation =
          new RandomAlignmentCreation(problem, 7, 0.2, new SeededRandom(1));

      // Act
      List<MsaSolution> solutions = creation.create();

      // Assert: 10 residues * 1.2 = 12 columns; toMsa() checks every row spells its protein
      assertEquals(7, solutions.size());
      assertEquals(12, creation.alignmentLength());
      for (MsaSolution solution : solutions) {
        assertEquals(12, solution.alignmentLength());
        assertEquals(2, solution.toMsa().size());
      }
    }

    @Test
    @DisplayName("given the same random seed, when creating twice, then the alignments are equal")
    void givenSameRandomSeed_whenCreatingTwice_thenAlignmentsAreEqual() {
      // Act
      List<MsaSolution> first =
          new RandomAlignmentCreation(problem, 3, 0.5, new SeededRandom(9)).create();
      List<MsaSolution> second =
          new RandomAlignmentCreation(problem, 3, 0.5, new SeededRandom(9)).create();
      List<MsaSolution> other =
          new RandomAlignmentCreation(problem, 3, 0.5, new SeededRandom(10)).create();

      // Assert
      assertEquals(Alignments.rows(first.get(2)), Alignments.rows(second.get(2)));
      assertNotEquals(
          first.stream().map(Alignments::rows).toList(),
          other.stream().map(Alignments::rows).toList());
    }
  }

  @Nested
  @DisplayName("When the arguments are invalid")
  class InvalidArguments {

    @ParameterizedTest
    @ValueSource(ints = {0, -1})
    @DisplayName("given a non-positive batch size, when creating the generator, then it fails")
    void givenNonPositiveBatchSize_whenCreatingGenerator_thenIllegalArgumentExceptionIsThrown(
        int batchSize) {
      // Act & Assert
      assertThrows(
          IllegalArgumentException.class,
          () -> new RandomAlignmentCreation(problem, batchSize, 0.2, new SeededRandom(1)));
    }

    @ParameterizedTest
    @ValueSource(doubles = {-0.1, Double.NaN, Double.POSITIVE_INFINITY})
    @DisplayName("given an invalid length increase, when creating the generator, then it fails")
    void givenInvalidLengthIncrease_whenCreatingGenerator_thenIllegalArgumentExceptionIsThrown(
        double increase) {
      // Act & Assert
      assertThrows(
          IllegalArgumentException.class,
          () -> new RandomAlignmentCreation(problem, 1, increase, new SeededRandom(1)));
    }
  }
}
