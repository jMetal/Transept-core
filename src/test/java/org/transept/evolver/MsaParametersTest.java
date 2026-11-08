package org.transept.evolver;

import static org.junit.jupiter.api.Assertions.assertInstanceOf;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.util.List;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.transept.model.ProteinSet;
import org.transept.operator.MergeAdjacentGapRunsMutation;
import org.transept.operator.RandomChoiceMutation;
import org.transept.operator.TmSafeBlockShiftMutation;
import org.transept.operator.TmSafeGapShiftMutation;
import org.transept.operator.TmSafeSinglePointCrossover;
import org.transept.operator.TmSafeUniformCrossover;
import org.transept.operator.TmSafeWindowRealignMutation;
import org.transept.testing.Alignments;
import org.transept.testing.SeededRandom;
import org.uma.evolver.parameter.type.CategoricalParameter;

@DisplayName("Unit tests for the operators and parameters of the MSA parameter spaces")
class MsaParametersTest {

  @ParameterizedTest
  @CsvSource({"tmSafeGapShift, true", "mergeGapRuns, false", "randomChoice, false"})
  @DisplayName("given a mutation name, when built, then the operator of that name is returned")
  void givenMutationName_whenBuilt_thenOperatorOfThatNameIsReturned(String name, boolean shift) {
    // Arrange
    ProteinSet proteins = Alignments.proteins("ACDEF:OOMMM", "ACDE:OOMM");

    // Act
    var mutation =
        MsaOperators.mutation(name, 0.5, 0, Alignments.unstructured(proteins), new SeededRandom(1));

    // Assert
    Class<?> expected =
        shift
            ? TmSafeGapShiftMutation.class
            : name.equals("mergeGapRuns")
                ? MergeAdjacentGapRunsMutation.class
                : RandomChoiceMutation.class;
    assertInstanceOf(expected, mutation);
  }

  @Test
  @DisplayName("given a limit, when building the mutation that chooses, then it is built")
  void givenLimit_whenBuildingRandomChoice_thenItIsBuilt() {
    // Arrange
    ProteinSet proteins = Alignments.proteins("ACDEF:OOMMM", "ACDE:OOMM");

    // Act & Assert
    assertInstanceOf(
        RandomChoiceMutation.class,
        MsaOperators.mutation(
            "randomChoice", 0.5, 4, Alignments.unstructured(proteins), new SeededRandom(1)));
  }

  @Test
  @DisplayName("given the new mutation names, when built, then their operators are returned")
  void givenNewMutationNames_whenBuilt_thenTheirOperatorsAreReturned() {
    // Arrange
    ProteinSet proteins = Alignments.proteins("ACDEF:OOMMM", "ACDE:OOMM");
    var structures = Alignments.unstructured(proteins);

    // Act & Assert
    assertInstanceOf(
        TmSafeBlockShiftMutation.class,
        MsaOperators.mutation("tmSafeBlockShift", 0.5, 0, structures, new SeededRandom(1)));
    assertInstanceOf(
        TmSafeWindowRealignMutation.class,
        MsaOperators.mutation("tmSafeWindowRealign", 0.5, 0, structures, new SeededRandom(1)));
  }

  @Test
  @DisplayName("given the crossover name, when built, then the TM-safe crossover is returned")
  void givenCrossoverName_whenBuilt_thenTmSafeCrossoverIsReturned() {
    // Act & Assert
    assertInstanceOf(
        TmSafeSinglePointCrossover.class,
        MsaOperators.crossover("tmSafeSinglePoint", 0.9, 1, new SeededRandom(1)));
    assertInstanceOf(
        TmSafeUniformCrossover.class,
        MsaOperators.crossover("tmSafeUniform", 0.9, 1, new SeededRandom(1)));
  }

  @Test
  @DisplayName("given unknown operator names, when built, then an exception is thrown")
  void givenUnknownOperatorNames_whenBuilt_thenIllegalArgumentExceptionIsThrown() {
    // Arrange
    ProteinSet proteins = Alignments.proteins("ACDEF:OOMMM", "ACDE:OOMM");

    // Act & Assert
    assertThrows(
        IllegalArgumentException.class,
        () -> MsaOperators.crossover("PMX", 0.9, 1, new SeededRandom(1)));
    assertThrows(
        IllegalArgumentException.class,
        () -> MsaOperators.mutation("swap", 0.1, 0, Alignments.unstructured(proteins), null));
  }

  @Test
  @DisplayName("given unknown values in a space, when creating the parameters, then they fail")
  void givenUnknownValuesInSpace_whenCreatingParameters_thenIllegalArgumentExceptionIsThrown() {
    // Arrange
    var factory = new MsaParameterFactory();

    // Act & Assert
    assertThrows(
        IllegalArgumentException.class, () -> factory.createParameter("crossover", List.of("PMX")));
    assertThrows(
        IllegalArgumentException.class, () -> factory.createParameter("mutation", List.of("swap")));
    assertThrows(
        IllegalArgumentException.class,
        () -> factory.createParameter("variation", List.of("differentialEvolution")));
    assertThrows(
        IllegalArgumentException.class,
        () -> factory.createParameter("createInitialSolutions", List.of("random")));
    assertThrows(IllegalArgumentException.class, () -> factory.createParameter("x", List.of()));
    assertInstanceOf(
        CategoricalParameter.class, factory.createParameter("algorithmResult", List.of("a")));
  }
}
