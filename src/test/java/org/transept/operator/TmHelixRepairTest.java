package org.transept.operator;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.List;
import java.util.Random;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.transept.alignment.Msa;
import org.transept.alignment.MsaSolution;
import org.transept.model.ProteinSet;
import org.transept.testing.Alignments;

@DisplayName("Unit tests for class TmHelixRepair")
class TmHelixRepairTest {

  private ProteinSet proteins;

  @BeforeEach
  void setUp() {
    // Helices: residues 2-5 of the first protein (CDEF), 2-4 of the second and third
    proteins =
        Alignments.proteins("ACDEFGHIK:OOMMMMOOO", "ACDEGHIK:OOMMMOOO", "ACEFGHIK:OOMMMOOO");
  }

  private static int breakingRuns(Msa msa) {
    int breaking = 0;
    for (int sequence = 0; sequence < msa.size(); sequence++) {
      breaking +=
          Alignments.helixBreakingRuns(msa.proteins().get(sequence), msa.rows().get(sequence));
    }
    return breaking;
  }

  @Test
  @DisplayName("given a gap inside a helix, when repairing, then it moves to the nearest boundary")
  void givenGapInsideHelix_whenRepairing_thenItMovesToTheNearestBoundary() {
    // Arrange: the gap of the second row lies after ACD, inside its helix DEG
    Msa msa =
        Msa.fromAlignedRows(proteins, List.of("ACDEFGHIK-", "ACD-EGHIK-", "ACEFGHIK--"));

    // Act
    TmHelixRepair.Result result = TmHelixRepair.repair(msa);

    // Assert: one residue back, just before the helix (one residue forward is still inside it)
    assertEquals(1, result.movedRuns());
    assertEquals("AC-DEGHIK-", result.alignment().alignedRow(1));
    assertEquals(msa.alignedRow(0), result.alignment().alignedRow(0));
    assertEquals(0, breakingRuns(result.alignment()));
    assertEquals(msa.length(), result.alignment().length());
  }

  @Test
  @DisplayName("given gaps only in loops and termini, when repairing, then nothing changes")
  void givenGapsOnlyInLoopsAndTermini_whenRepairing_thenNothingChanges() {
    // Arrange
    Msa msa =
        Msa.fromAlignedRows(proteins, List.of("-ACDEFGHIK", "-AC-DEGHIK", "--ACEFGHIK"));

    // Act
    TmHelixRepair.Result result = TmHelixRepair.repair(msa);

    // Assert
    assertEquals(0, result.movedRuns());
    assertEquals(msa.rows(), result.alignment().rows());
  }

  @ParameterizedTest
  @ValueSource(ints = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10})
  @DisplayName("given random alignments, when repairing, then no helix is broken and rows keep")
  void givenRandomAlignments_whenRepairing_thenNoHelixIsBrokenAndRowsKeepTheirResidues(int seed) {
    // Arrange
    MsaSolution solution = Alignments.randomSolution(proteins, 6, new Random(seed));
    Msa msa = solution.toMsa();

    // Act
    Msa repaired = TmHelixRepair.repair(msa).alignment();

    // Assert
    assertEquals(0, breakingRuns(repaired));
    assertEquals(msa.length(), repaired.length());
    for (int sequence = 0; sequence < msa.size(); sequence++) {
      assertEquals(
          msa.alignedRow(sequence).replace("-", ""),
          repaired.alignedRow(sequence).replace("-", ""));
    }
  }
}
