package org.tmm2structalign.algorithm;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.times;
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
import org.tmm2structalign.objective.TopologyAwareSumOfPairs;
import org.tmm2structalign.operator.TmSafeGapShiftMutation;
import org.tmm2structalign.problem.StructuralMsaProblem;
import org.tmm2structalign.testing.Alignments;
import org.tmm2structalign.testing.SeededRandom;
import org.uma.jmetal.operator.mutation.MutationOperator;

@DisplayName("Unit tests for class PerturbedSeedCreation")
class PerturbedSeedCreationTest {

  private ProteinSet proteins;
  private StructuralMsaProblem problem;
  private List<Msa> seeds;

  @BeforeEach
  void setUp() {
    proteins = Alignments.proteins("ACDEFGHIK:OOMMMMOOO", "ACDEGHIK:OOMMMOOO");
    problem =
        new StructuralMsaProblem(
            "p", proteins, List.of(TopologyAwareSumOfPairs.withDefaultParameters()));
    seeds =
        List.of(
            Msa.fromAlignedRows(proteins, List.of("ACDEFGHIK--", "ACDE-GHIK--")),
            Msa.fromAlignedRows(proteins, List.of("--ACDEFGHIK", "--ACDEGHIK-")));
  }

  @Nested
  @DisplayName("When perturbing seeds")
  class Perturbation {

    @Test
    @DisplayName("given no mutations, when creating, then every solution is one of the seeds")
    void givenNoMutations_whenCreating_thenEverySolutionIsOneOfTheSeeds() {
      // Arrange
      @SuppressWarnings("unchecked")
      MutationOperator<MsaSolution> mutation = mock(MutationOperator.class);
      PerturbedSeedCreation creation =
          new PerturbedSeedCreation(problem, seeds, 6, 0, mutation, new SeededRandom(1));

      // Act
      List<MsaSolution> solutions = creation.create();

      // Assert
      assertEquals(6, solutions.size());
      for (MsaSolution solution : solutions) {
        assertTrue(
            seeds.stream().anyMatch(seed -> seed.rows().equals(solution.variables())));
      }
    }

    @Test
    @DisplayName("given k mutations, when creating, then each solution is mutated k times")
    void givenMutations_whenCreating_thenEachSolutionIsMutatedThatManyTimes() {
      // Arrange
      @SuppressWarnings("unchecked")
      MutationOperator<MsaSolution> mutation = mock(MutationOperator.class);
      when(mutation.execute(any())).thenAnswer(invocation -> invocation.getArgument(0));
      PerturbedSeedCreation creation =
          new PerturbedSeedCreation(problem, seeds, 4, 3, mutation, new SeededRandom(1));

      // Act
      creation.create();

      // Assert
      verify(mutation, times(12)).execute(any());
    }

    @Test
    @DisplayName("given the real mutation, when creating, then alignments stay valid")
    void givenRealMutation_whenCreating_thenAlignmentsStayValid() {
      // Arrange
      PerturbedSeedCreation creation =
          new PerturbedSeedCreation(
              problem,
              seeds,
              20,
              PerturbedSeedCreation.DEFAULT_MUTATIONS,
              new TmSafeGapShiftMutation(
                  1.0, Alignments.fullyStructured(proteins), 60, new SeededRandom(2)),
              new SeededRandom(3));

      // Act
      List<MsaSolution> solutions = creation.create();

      // Assert
      for (MsaSolution solution : solutions) {
        assertEquals(11, solution.alignmentLength());
        assertEquals(2, solution.toMsa().size());
      }
    }
  }

  @Nested
  @DisplayName("When the arguments are invalid")
  class InvalidArguments {

    @Test
    @DisplayName("given invalid arguments, when creating the generator, then it fails")
    void givenInvalidArguments_whenCreatingGenerator_thenIllegalArgumentExceptionIsThrown() {
      // Arrange
      @SuppressWarnings("unchecked")
      MutationOperator<MsaSolution> mutation = mock(MutationOperator.class);
      SeededRandom random = new SeededRandom(1);

      // Act & Assert
      assertThrows(
          IllegalArgumentException.class,
          () -> new PerturbedSeedCreation(problem, List.of(), 1, 1, mutation, random));
      assertThrows(
          IllegalArgumentException.class,
          () -> new PerturbedSeedCreation(problem, seeds, 0, 1, mutation, random));
      assertThrows(
          IllegalArgumentException.class,
          () -> new PerturbedSeedCreation(problem, seeds, 1, -1, mutation, random));
    }
  }
}
