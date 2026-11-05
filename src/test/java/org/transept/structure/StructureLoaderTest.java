package org.transept.structure;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Optional;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.transept.model.ProteinSequence;
import org.transept.model.ProteinSet;
import org.transept.testing.TestFiles;

@DisplayName("Unit tests for classes StructureLoader and ProteinStructures")
class StructureLoaderTest {

  @TempDir Path directory;

  private StructureLoader loader;
  private ProteinSet proteins;

  @BeforeEach
  void setUp() throws IOException {
    loader = new StructureLoader(MappingPolicy.VERIFIED);
    proteins =
        new ProteinSet(
            List.of(
                ProteinSequence.withoutTopology("good", "ACDEF"),
                ProteinSequence.withoutTopology("missing", "ACDEF"),
                ProteinSequence.withoutTopology("broken", "ACDEF"),
                ProteinSequence.withoutTopology("extended", "MMACDEF")));
    TestFiles.writeLinearStructure(directory, "good", 5, 3.8f);
    TestFiles.writePdb(directory.resolve("good.pdb"), "ACDEF");
    Files.writeString(directory.resolve("broken_D.csv"), "0,1\n2,0\n");
    Files.writeString(directory.resolve("broken_idx.csv"), "1,2");
    TestFiles.writeLinearStructure(directory, "extended", 5, 3.8f);
    TestFiles.writePdb(directory.resolve("extended.pdb"), "ACDEF");
  }

  @Nested
  @DisplayName("When loading a dataset")
  class Loading {

    @Test
    @DisplayName("given mixed structural files, when loading, then every protein gets a status")
    void givenMixedFiles_whenLoading_thenEveryProteinGetsStatus() throws IOException {
      // Act
      ProteinStructures structures = loader.load(proteins, directory, Optional.of(directory));

      // Assert
      assertEquals(4, structures.size());
      assertEquals(MappingStatus.NUMBERED, structures.mapping(0).status());
      assertEquals(MappingStatus.MISSING_FILES, structures.mapping(1).status());
      assertEquals(MappingStatus.INVALID_FILES, structures.mapping(2).status());
      assertEquals(MappingStatus.ALIGNED, structures.mapping(3).status());
      assertEquals(List.of("missing", "broken"), structures.proteinsWithoutStructure());
      assertEquals(4, structures.mappings().size());
    }

    @Test
    @DisplayName("given a usable structure, when querying, then distances and rows are exposed")
    void givenUsableStructure_whenQuerying_thenDistancesAndRowsAreExposed() throws IOException {
      // Act
      ProteinStructures structures = loader.load(proteins, directory, Optional.of(directory));

      // Assert
      assertEquals(3.8f, structures.distances(0).orElseThrow()[0][1], 1e-6f);
      assertTrue(structures.distances(1).isEmpty());
      assertEquals(0, structures.structureIndex(3, 2));
      assertEquals(-1, structures.structureIndex(3, 0));
      assertEquals(-1, structures.structureIndex(1, 0));
    }

    @Test
    @DisplayName("given no PDB directory, when loading, then mappings are unverified")
    void givenNoPdbDirectory_whenLoading_thenMappingsAreUnverified() throws IOException {
      // Act
      ProteinStructures structures = loader.load(proteins, directory, Optional.empty());

      // Assert
      assertEquals(MappingStatus.UNVERIFIED, structures.mapping(0).status());
      assertEquals(MappingStatus.UNVERIFIED, structures.mapping(3).status());
    }
  }

  @Nested
  @DisplayName("When building the collection directly")
  class Construction {

    @Test
    @DisplayName("given lists of different sizes, when creating, then an exception is thrown")
    void givenDifferentSizes_whenCreating_thenIllegalArgumentExceptionIsThrown() {
      // Arrange
      List<StructureMapping> mappings =
          List.of(StructureMapping.unusable("p", MappingStatus.MISSING_FILES, ""));

      // Act & Assert
      assertThrows(
          IllegalArgumentException.class, () -> new ProteinStructures(mappings, List.of()));
    }
  }
}
