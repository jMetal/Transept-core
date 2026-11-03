package org.tmm2structalign.structure;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Optional;
import java.util.stream.IntStream;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.tmm2structalign.model.ProteinSequence;

@DisplayName("Unit tests for class SequenceStructureMapper")
class SequenceStructureMapperTest {

  private static final String CORE = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKL";

  @Nested
  @DisplayName("When using the verified policy")
  class Verified {

    private final SequenceStructureMapper mapper =
        new SequenceStructureMapper(MappingPolicy.VERIFIED);

    @Test
    @DisplayName("given identical sequences, when mapping, then the numbering is used")
    void givenIdenticalSequences_whenMapping_thenNumberingIsUsed() {
      // Arrange
      ProteinSequence protein = ProteinSequence.withoutTopology("p", CORE);

      // Act
      StructureMapping mapping = mapper.map(protein, numbers(CORE.length()), Optional.of(CORE));

      // Assert
      assertEquals(MappingStatus.NUMBERED, mapping.status());
      assertEquals(1.0, mapping.identity());
      assertEquals(1.0, mapping.coverage());
      assertArrayEquals(IntStream.range(0, CORE.length()).toArray(), mapping.residueToStructure());
      assertTrue(mapping.explanation().isEmpty());
    }

    @Test
    @DisplayName("given a point variant, when mapping, then the numbering is still accepted")
    void givenPointVariant_whenMapping_thenNumberingIsAccepted() {
      // Arrange
      String variant = "W" + CORE.substring(1);
      ProteinSequence protein = ProteinSequence.withoutTopology("p", CORE);

      // Act
      StructureMapping mapping =
          mapper.map(protein, numbers(CORE.length()), Optional.of(variant));

      // Assert
      assertEquals(MappingStatus.NUMBERED, mapping.status());
      assertTrue(mapping.identity() < 1.0);
    }

    @Test
    @DisplayName("given an N-terminal extension (OXER1 case), when mapping, then it aligns")
    void givenAminoTerminalExtension_whenMapping_thenAlignmentIsUsed() {
      // Arrange
      String extension = "MSTPQ";
      ProteinSequence protein = ProteinSequence.withoutTopology("oxer1", extension + CORE);

      // Act
      StructureMapping mapping = mapper.map(protein, numbers(CORE.length()), Optional.of(CORE));

      // Assert
      assertEquals(MappingStatus.ALIGNED, mapping.status());
      assertEquals(1.0, mapping.identity());
      assertEquals(-1, mapping.residueToStructure()[0]);
      assertEquals(0, mapping.residueToStructure()[extension.length()]);
      assertEquals(CORE.length() / (double) protein.length(), mapping.coverage(), 1e-12);
      assertTrue(mapping.explanation().orElseThrow().contains("alignment identity"));
    }

    @Test
    @DisplayName("given an unrelated structure, when mapping, then the structure is rejected")
    void givenUnrelatedStructure_whenMapping_thenStructureIsRejected() {
      // Arrange
      ProteinSequence protein = ProteinSequence.withoutTopology("p", "AAAAAAAAAA");

      // Act
      StructureMapping mapping = mapper.map(protein, numbers(10), Optional.of("WWWWWWWWWW"));

      // Assert
      assertEquals(MappingStatus.LOW_IDENTITY, mapping.status());
      assertEquals(0, mapping.residueToStructure().length);
    }

    @Test
    @DisplayName("given residue numbers beyond the sequence, when mapping, then it aligns")
    void givenNumbersBeyondSequence_whenMapping_thenAlignmentIsUsed() {
      // Arrange
      ProteinSequence protein = ProteinSequence.withoutTopology("p", CORE);
      int[] shifted = IntStream.rangeClosed(11, 10 + CORE.length()).toArray();

      // Act
      StructureMapping mapping = mapper.map(protein, shifted, Optional.of(CORE));

      // Assert
      assertEquals(MappingStatus.ALIGNED, mapping.status());
      assertEquals(0, mapping.residueToStructure()[0]);
    }

    @Test
    @DisplayName("given no PDB sequence, when mapping, then the numbering is unverified")
    void givenNoPdbSequence_whenMapping_thenNumberingIsUnverified() {
      // Arrange
      ProteinSequence protein = ProteinSequence.withoutTopology("p", CORE);

      // Act
      StructureMapping mapping = mapper.map(protein, numbers(5), Optional.empty());

      // Assert
      assertEquals(MappingStatus.UNVERIFIED, mapping.status());
      assertEquals(5 / (double) CORE.length(), mapping.coverage(), 1e-12);
    }

    @Test
    @DisplayName("given a PDB sequence of the wrong length, when mapping, then it is rejected")
    void givenPdbSequenceOfWrongLength_whenMapping_thenIllegalArgumentExceptionIsThrown() {
      // Arrange
      ProteinSequence protein = ProteinSequence.withoutTopology("p", CORE);

      // Act & Assert
      assertThrows(
          IllegalArgumentException.class,
          () -> mapper.map(protein, numbers(3), Optional.of("AC")));
    }
  }

  @Nested
  @DisplayName("When using the other policies")
  class OtherPolicies {

    @Test
    @DisplayName("given an extension and numbering-only policy, when mapping, then it is rejected")
    void givenExtensionAndNumberingOnly_whenMapping_thenStructureIsRejected() {
      // Arrange
      SequenceStructureMapper mapper = new SequenceStructureMapper(MappingPolicy.NUMBERED_ONLY);
      ProteinSequence protein = ProteinSequence.withoutTopology("p", "MSTPQ" + CORE);

      // Act
      StructureMapping mapping = mapper.map(protein, numbers(CORE.length()), Optional.of(CORE));

      // Assert
      assertEquals(MappingStatus.LOW_IDENTITY, mapping.status());
    }

    @Test
    @DisplayName("given an extension and the legacy policy, when mapping, then it is mis-mapped")
    void givenExtensionAndLegacyPolicy_whenMapping_thenOriginalMisMappingIsReproduced() {
      // Arrange
      SequenceStructureMapper mapper = new SequenceStructureMapper(MappingPolicy.LEGACY);
      ProteinSequence protein = ProteinSequence.withoutTopology("p", "MSTPQ" + CORE);

      // Act
      StructureMapping mapping = mapper.map(protein, numbers(CORE.length()), Optional.of(CORE));

      // Assert
      assertEquals(MappingStatus.UNVERIFIED, mapping.status());
      assertEquals(0, mapping.residueToStructure()[0]);
    }

    @Test
    @DisplayName("given numbers beyond the sequence and the legacy policy, then files are invalid")
    void givenNumbersBeyondSequenceAndLegacyPolicy_whenMapping_thenFilesAreInvalid() {
      // Arrange
      SequenceStructureMapper mapper = new SequenceStructureMapper(MappingPolicy.LEGACY);
      ProteinSequence protein = ProteinSequence.withoutTopology("p", "ACD");

      // Act
      StructureMapping mapping = mapper.map(protein, numbers(4), Optional.empty());

      // Assert
      assertEquals(MappingStatus.INVALID_FILES, mapping.status());
    }

    @Test
    @DisplayName("given repeated residue numbers, when mapping, then the first row is kept")
    void givenRepeatedResidueNumbers_whenMapping_thenFirstRowIsKept() {
      // Arrange
      SequenceStructureMapper mapper = new SequenceStructureMapper(MappingPolicy.LEGACY);
      ProteinSequence protein = ProteinSequence.withoutTopology("p", "ACD");

      // Act
      StructureMapping mapping = mapper.map(protein, new int[] {1, 2, 2, 3}, Optional.empty());

      // Assert
      assertArrayEquals(new int[] {0, 1, 3}, mapping.residueToStructure());
    }
  }

  @Nested
  @DisplayName("When aligning sequences")
  class Alignment {

    @Test
    @DisplayName("given an internal insertion, when aligning, then residues around it map")
    void givenInternalInsertion_whenAligning_thenResiduesAroundItMap() {
      // Act
      int[] mapping = SequenceStructureMapper.semiGlobalAlignment("ACDEXXXFGHIK", "ACDEFGHIK");

      // Assert
      assertEquals(3, mapping[3]);
      assertEquals(4, mapping[7]);
      assertEquals(8, mapping[11]);
    }
  }

  private static int[] numbers(int count) {
    return IntStream.rangeClosed(1, count).toArray();
  }
}
