package org.tmm2structalign.problem;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.verify;
import static org.mockito.Mockito.when;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.tmm2structalign.algorithm.AlgorithmConfiguration;
import org.tmm2structalign.alignment.Msa;
import org.tmm2structalign.alignment.MsaSolution;
import org.tmm2structalign.model.ProteinSet;
import org.tmm2structalign.objective.Objective;
import org.tmm2structalign.run.Experiment;
import org.tmm2structalign.run.RunConfiguration;
import org.tmm2structalign.structure.MappingPolicy;
import org.tmm2structalign.testing.Alignments;
import org.tmm2structalign.testing.SyntheticDataset;

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

  @Nested
  @DisplayName("When evaluating concurrently")
  class Concurrency {

    @TempDir Path directory;

    @Test
    @DisplayName("given many threads, when evaluating copies, then results match sequential ones")
    void givenManyThreads_whenEvaluatingCopies_thenResultsEqualSequentialOnes() throws Exception {
      // Arrange: a real problem (SOPwTP + lDDT with structures) shared by every thread
      Path root = SyntheticDataset.write(directory.resolve("GPCRdb"));
      StructuralMsaProblem shared =
          Experiment.load(
                  new RunConfiguration(
                      root,
                      SyntheticDataset.NAME,
                      SyntheticDataset.SEED_METHODS,
                      directory.resolve("unused"),
                      1L,
                      AlgorithmConfiguration.withDefaults(10, 4, 1),
                      MappingPolicy.VERIFIED,
                      15f,
                      10))
              .problem();
      Random random = new Random(17);
      List<MsaSolution> originals = new ArrayList<>();
      for (int index = 0; index < 64; index++) {
        originals.add(Alignments.randomSolution(shared.proteins(), 40, random));
      }
      List<double[]> expected = new ArrayList<>();
      for (MsaSolution original : originals) {
        expected.add(shared.evaluate(original.copy()).objectives().clone());
      }
      ExecutorService pool = Executors.newFixedThreadPool(8);

      try {
        for (int round = 0; round < 5; round++) {
          // Act
          List<MsaSolution> copies = originals.stream().map(MsaSolution::copy).toList();
          List<Future<MsaSolution>> futures = new ArrayList<>();
          for (MsaSolution copy : copies) {
            futures.add(pool.submit(() -> shared.evaluate(copy)));
          }

          // Assert
          for (int index = 0; index < futures.size(); index++) {
            assertArrayEquals(expected.get(index), futures.get(index).get().objectives());
          }
        }
      } finally {
        pool.shutdownNow();
      }
    }
  }
}
