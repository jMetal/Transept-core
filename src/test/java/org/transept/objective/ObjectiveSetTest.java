package org.transept.objective;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.List;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.transept.model.ProteinSet;
import org.transept.testing.Alignments;

@DisplayName("Unit tests for enum ObjectiveSet")
class ObjectiveSetTest {

  @ParameterizedTest
  @CsvSource({
    "LDDT, lDDT",
    "SOPWTP_LDDT, SOPwTP lDDT",
    "SOPWTP_LDDT_COVERAGE, SOPwTP lDDT coverage"
  })
  @DisplayName("given a formulation, when creating it, then its objectives come in output order")
  void givenFormulation_whenCreating_thenObjectivesComeInOutputOrder(
      ObjectiveSet set, String names) {
    // Arrange
    ProteinSet proteins = Alignments.proteins("ACDEF:OOMMM", "ACDE:OOMM");

    // Act
    List<Objective> objectives =
        set.create(proteins, Alignments.unstructured(proteins), StructuralLddt.DEFAULT_CUTOFF);

    // Assert
    assertEquals(List.of(names.split(" ")), objectives.stream().map(Objective::name).toList());
  }
}
