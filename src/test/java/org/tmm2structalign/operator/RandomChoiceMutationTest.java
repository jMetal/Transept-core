package org.tmm2structalign.operator;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.never;
import static org.mockito.Mockito.times;
import static org.mockito.Mockito.verify;
import static org.mockito.Mockito.when;

import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.tmm2structalign.alignment.MsaSolution;
import org.tmm2structalign.model.ProteinSet;
import org.tmm2structalign.testing.Alignments;
import org.uma.jmetal.operator.mutation.MutationOperator;

@DisplayName("Unit tests for class RandomChoiceMutation")
class RandomChoiceMutationTest {

  private RandomSource random;
  private MutationOperator<MsaSolution> first;
  private MutationOperator<MsaSolution> second;
  private MsaSolution solution;

  @BeforeEach
  @SuppressWarnings("unchecked")
  void setUp() {
    random = mock(RandomSource.class);
    first = mock(MutationOperator.class);
    second = mock(MutationOperator.class);
    ProteinSet proteins = Alignments.proteins("ACD:OOO", "ACD:OOO");
    solution = Alignments.solution(proteins, "AC-D", "ACD-");
    when(first.execute(any())).thenAnswer(invocation -> invocation.getArgument(0));
    when(second.execute(any())).thenAnswer(invocation -> invocation.getArgument(0));
  }

  @Test
  @DisplayName("given three mutations per application, when mutating, then the drawn ones run")
  void givenThreeMutationsPerApplication_whenMutating_thenTheDrawnOnesRun() {
    // Arrange
    when(random.nextDouble()).thenReturn(0.0);
    when(random.nextInt(0, 1)).thenReturn(1, 0, 1);
    RandomChoiceMutation mutation =
        new RandomChoiceMutation(1.0, List.of(first, second), 3, random);

    // Act
    MsaSolution result = mutation.execute(solution);

    // Assert
    assertSame(solution, result);
    verify(first, times(1)).execute(solution);
    verify(second, times(2)).execute(solution);
  }

  @Test
  @DisplayName("given the probability is not reached, when mutating, then no operator runs")
  void givenProbabilityIsNotReached_whenMutating_thenNoOperatorRuns() {
    // Arrange
    when(random.nextDouble()).thenReturn(0.9);
    RandomChoiceMutation mutation =
        new RandomChoiceMutation(0.5, List.of(first, second), 1, random);

    // Act
    mutation.execute(solution);

    // Assert
    verify(first, never()).execute(any());
    verify(second, never()).execute(any());
    assertEquals(0.5, mutation.mutationProbability());
    assertEquals(List.of(first, second), mutation.operators());
  }

  @Test
  @DisplayName("given invalid arguments, when creating or mutating, then it fails")
  void givenInvalidArguments_whenCreatingOrMutating_thenIllegalArgumentExceptionIsThrown() {
    assertThrows(
        IllegalArgumentException.class,
        () -> new RandomChoiceMutation(1.0, List.of(), 1, random));
    assertThrows(
        IllegalArgumentException.class,
        () -> new RandomChoiceMutation(1.0, List.of(first), 0, random));
    assertThrows(
        IllegalArgumentException.class,
        () -> new RandomChoiceMutation(2.0, List.of(first), 1, random));
    RandomChoiceMutation mutation = new RandomChoiceMutation(1.0, List.of(first), 1, random);
    assertThrows(IllegalArgumentException.class, () -> mutation.execute(null));
  }
}
