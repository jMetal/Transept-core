package org.tmm2structalign.problem;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.verify;
import static org.mockito.Mockito.when;

import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.tmm2structalign.alignment.Msa;
import org.tmm2structalign.alignment.MsaSolution;
import org.tmm2structalign.model.ProteinSet;
import org.tmm2structalign.objective.Objective;
import org.tmm2structalign.testing.Alignments;

@DisplayName("Unit tests for class StructuralMsaProblem")
class StructuralMsaProblemTest {

  private ProteinSet proteins;
  private Objective maximized;
  private Objective minimized;
  private StructuralMsaProblem problem;

  @BeforeEach
  void setUp() {
    proteins = Alignments.proteins("ACD:OOO", "AD:OO");
    maximized = mock(Objective.class);
    minimized = mock(Objective.class);
    when(maximized.isMaximized()).thenReturn(true);
    when(minimized.isMaximized()).thenReturn(false);
    when(maximized.evaluate(any())).thenReturn(5.0);
    when(minimized.evaluate(any())).thenReturn(2.0);
    problem = new StructuralMsaProblem("dataset", proteins, List.of(maximized, minimized));
  }

  @Nested
  @DisplayName("When evaluating solutions")
  class Evaluation {

    @Test
    @DisplayName("given gap-only columns, when evaluating, then they are removed first")
    void givenGapOnlyColumns_whenEvaluating_thenTheyAreRemovedFirst() {
      // Arrange
      MsaSolution solution = Alignments.solution(proteins, "-ACD-", "-A-D-");

      // Act
      problem.evaluate(solution);

      // Assert
      assertEquals(3, solution.alignmentLength());
      verify(maximized).evaluate(any(Msa.class));
    }

    @Test
    @DisplayName("given objectives, when evaluating, then maximised values are negated")
    void givenObjectives_whenEvaluating_thenMaximisedValuesAreNegated() {
      // Arrange
      MsaSolution solution = Alignments.solution(proteins, "ACD", "A-D");

      // Act
      problem.evaluate(solution);

      // Assert
      assertEquals(-5.0, solution.objectives()[0]);
      assertEquals(2.0, solution.objectives()[1]);
      assertEquals(5.0, problem.naturalValue(solution, 0));
      assertEquals(2.0, problem.naturalValue(solution, 1));
    }
  }

  @Nested
  @DisplayName("When creating solutions")
  class Creation {

    @Test
    @DisplayName("given an alignment of the problem, when creating, then a solution is returned")
    void givenAlignmentOfProblem_whenCreating_thenSolutionIsReturned() {
      // Arrange
      Msa msa = Msa.fromAlignedRows(proteins, List.of("ACD", "A-D"));

      // Act
      MsaSolution solution = problem.createSolution(msa);

      // Assert
      assertEquals(msa.rows(), solution.variables());
      assertEquals(2, solution.objectives().length);
    }

    @Test
    @DisplayName("given an alignment of other proteins, when creating, then it is rejected")
    void givenAlignmentOfOtherProteins_whenCreating_thenIllegalArgumentExceptionIsThrown() {
      // Arrange
      ProteinSet others = Alignments.proteins("ACD:OOO", "AD:OO");
      Msa msa = Msa.fromAlignedRows(others, List.of("ACD", "A-D"));

      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> problem.createSolution(msa));
    }

    @Test
    @DisplayName("given no alignment, when creating a random solution, then it is unsupported")
    void givenNoAlignment_whenCreatingRandomSolution_thenUnsupportedOperationExceptionIsThrown() {
      // Act & Assert
      assertThrows(UnsupportedOperationException.class, () -> problem.createSolution());
    }
  }

  @Nested
  @DisplayName("When describing the problem")
  class Description {

    @Test
    @DisplayName("given a problem, when querying sizes, then they reflect proteins and objectives")
    void givenProblem_whenQueryingSizes_thenTheyReflectProteinsAndObjectives() {
      // Act & Assert
      assertEquals("dataset", problem.name());
      assertEquals(2, problem.numberOfVariables());
      assertEquals(2, problem.numberOfObjectives());
      assertEquals(0, problem.numberOfConstraints());
      assertEquals(proteins, problem.proteins());
      assertEquals(2, problem.objectives().size());
    }

    @Test
    @DisplayName("given no objectives, when creating the problem, then it is rejected")
    void givenNoObjectives_whenCreatingProblem_thenIllegalArgumentExceptionIsThrown() {
      // Act & Assert
      assertThrows(
          IllegalArgumentException.class,
          () -> new StructuralMsaProblem("dataset", proteins, List.of()));
    }
  }
}
