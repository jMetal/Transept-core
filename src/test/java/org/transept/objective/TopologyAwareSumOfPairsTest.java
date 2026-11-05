package org.transept.objective;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.transept.alignment.Msa;
import org.transept.model.ProteinSequence;
import org.transept.model.ProteinSet;

@DisplayName("Unit tests for classes TopologyAwareSumOfPairs, SubstitutionMatrix and GapPenalties")
class TopologyAwareSumOfPairsTest {

  private TopologyAwareSumOfPairs objective;

  @BeforeEach
  void setUp() {
    objective = TopologyAwareSumOfPairs.withDefaultParameters();
  }

  @Nested
  @DisplayName("When scoring hand-computed alignments")
  class HandComputed {

    @ParameterizedTest(name = "{0}")
    @CsvSource(
        delimiter = '|',
        value = {
          // description | residues/topology 1 | residues/topology 2 | aligned 1 | aligned 2 | score
          "non-TM matches use BLOSUM62 (4 + 9)        | AC:OO  | AC:OO   | AC    | AC    | 13",
          "TM-TM matches use PHAT (5 + 7)             | AC:MM  | AC:MM   | AC    | AC    | 12",
          "mixed TM/non-TM pairs use BLOSUM62         | AC:MM  | AC:OO   | AC    | AC    | 13",
          "affine gap: open 3 + extend 1              | AC:OO  | AGGC:OOOO | A--C  | AGGC  | 9",
          "gap-gap column skipped, gap stays open     | AC:OO  | AGGC:OOOO | A---C | AG-GC | 9",
          "gaps on alternating sides both open        | AC:OO  | GC:OO   | A-C   | -GC   | 3",
          "gap between TM residues uses TM open (8)   | AC:MM  | AGC:OOO | A-C   | AGC   | 5",
          "only gap-gap columns besides matches       | AC:OO  | AC:OO   | A-C   | A-C   | 13"
        })
    @DisplayName("given an alignment, when scoring, then the hand-computed value is returned")
    void givenAlignment_whenScoring_thenHandComputedValueIsReturned(
        String description,
        String first,
        String second,
        String alignedFirst,
        String alignedSecond,
        double expected) {
      // Arrange
      Msa msa = alignment(first, second, alignedFirst, alignedSecond);

      // Act
      double score = objective.evaluate(msa);

      // Assert
      assertEquals(expected, score, description);
    }

    private Msa alignment(
        String first, String second, String alignedFirst, String alignedSecond) {
      ProteinSet proteins = new ProteinSet(List.of(protein("a", first), protein("b", second)));
      return Msa.fromAlignedRows(proteins, List.of(alignedFirst.trim(), alignedSecond.trim()));
    }

    private ProteinSequence protein(String id, String residuesAndTopology) {
      String[] parts = residuesAndTopology.trim().split(":");
      return new ProteinSequence(id, parts[0], parts[1]);
    }
  }

  @Nested
  @DisplayName("When checking general properties")
  class Properties {

    @Test
    @DisplayName("given three sequences in any order, when scoring, then the score is the same")
    void givenSequencesInAnyOrder_whenScoring_thenScoreIsTheSame() {
      // Arrange
      ProteinSequence a = new ProteinSequence("a", "ACDEFG", "OMMMMO");
      ProteinSequence b = new ProteinSequence("b", "ACEFG", "OMMMO");
      ProteinSequence c = new ProteinSequence("c", "CDEG", "MMMO");
      Msa forward =
          Msa.fromAlignedRows(
              new ProteinSet(List.of(a, b, c)), List.of("ACDEFG", "AC-EFG", "-CDE-G"));
      Msa reversed =
          Msa.fromAlignedRows(
              new ProteinSet(List.of(c, b, a)), List.of("-CDE-G", "AC-EFG", "ACDEFG"));

      // Act & Assert
      assertEquals(objective.evaluate(forward), objective.evaluate(reversed));
      assertEquals("SOPwTP", objective.name());
      assertTrue(objective.isMaximized());
    }
  }

  @Nested
  @DisplayName("When using substitution matrices and gap penalties")
  class MatricesAndPenalties {

    @Test
    @DisplayName("given standard residues, when looking up, then matrices are symmetric")
    void givenStandardResidues_whenLookingUp_thenMatricesAreSymmetric() {
      // Arrange
      String residues = "ARNDCQEGHILKMFPSTWYV";

      // Act & Assert
      for (char first : residues.toCharArray()) {
        for (char second : residues.toCharArray()) {
          assertEquals(
              SubstitutionMatrix.PHAT.score(first, second),
              SubstitutionMatrix.PHAT.score(second, first));
          assertEquals(
              SubstitutionMatrix.BLOSUM62.score(first, second),
              SubstitutionMatrix.BLOSUM62.score(second, first));
        }
      }
      assertEquals(11, SubstitutionMatrix.BLOSUM62.score('W', 'W'));
      assertEquals(-10, SubstitutionMatrix.PHAT.score('C', 'K'));
      assertEquals("PHAT", SubstitutionMatrix.PHAT.name());
    }

    @Test
    @DisplayName("given ambiguity codes, when looking up BLOSUM62, then the typos are fixed")
    void givenAmbiguityCodes_whenLookingUp_thenOriginalTyposAreFixed() {
      // Act & Assert
      assertEquals(3, SubstitutionMatrix.BLOSUM62.score('B', 'N'));
      assertEquals(3, SubstitutionMatrix.BLOSUM62.score('Z', 'Q'));
    }

    @ParameterizedTest
    @CsvSource({"PHAT, X", "PHAT, B", "BLOSUM62, U", "BLOSUM62, a", "BLOSUM62, é"})
    @DisplayName("given a residue outside the alphabet, when looking up, then it is rejected")
    void givenResidueOutsideAlphabet_whenLookingUp_thenIllegalArgumentExceptionIsThrown(
        String matrixName, char residue) {
      // Arrange
      SubstitutionMatrix matrix =
          "PHAT".equals(matrixName) ? SubstitutionMatrix.PHAT : SubstitutionMatrix.BLOSUM62;

      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> matrix.score(residue, 'A'));
    }

    @ParameterizedTest
    @CsvSource({"-1, 1", "1, -1", "NaN, 1", "1, Infinity"})
    @DisplayName("given invalid penalties, when creating them, then an exception is thrown")
    void givenInvalidPenalties_whenCreating_thenIllegalArgumentExceptionIsThrown(
        double open, double extend) {
      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> new GapPenalties(open, extend));
    }
  }
}
