package org.tmm2structalign.io;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.tmm2structalign.model.ProteinSet;
import org.tmm2structalign.model.Topology;

@DisplayName("Unit tests for class ThreeLineTopologyReader")
class ThreeLineTopologyReaderTest {

  @TempDir Path directory;

  private ThreeLineTopologyReader reader;

  @BeforeEach
  void setUp() {
    reader = new ThreeLineTopologyReader();
  }

  @Nested
  @DisplayName("When reading well-formed files")
  class WellFormed {

    @Test
    @DisplayName("given annotated headers, when reading, then ids, residues and topology match")
    void givenAnnotatedHeaders_whenReading_thenFieldsMatch() throws IOException {
      // Arrange
      Path file =
          write(">fzd2_human | SP+TM\nacdef\nOMMMI\n\n>smo_human\nGHIK\nSSOU\n");

      // Act
      ProteinSet proteins = reader.read(file);

      // Assert
      assertEquals(2, proteins.size());
      assertEquals("fzd2_human", proteins.get(0).id());
      assertEquals("ACDEF", proteins.get(0).residues());
      assertEquals(Topology.TRANSMEMBRANE, proteins.get(0).topologyAt(1));
      assertEquals("smo_human", proteins.get(1).id());
      assertEquals(Topology.SIGNAL_PEPTIDE, proteins.get(1).topologyAt(0));
    }
  }

  @Nested
  @DisplayName("When reading malformed files")
  class Malformed {

    @ParameterizedTest
    @ValueSource(
        strings = {
          ">a\nAC\nOO\n>b\nAC\n", // incomplete record
          "a\nAC\nOO\n>b\nAC\nOO\n", // header without '>'
          ">a\nACD\nOO\n>b\nAC\nOO\n", // topology shorter than sequence
          ">a\nAC\nOO\n>a\nAC\nOO\n" // duplicated id
        })
    @DisplayName("given a malformed file, when reading, then a descriptive exception is thrown")
    void givenMalformedFile_whenReading_thenInvalidInputFileExceptionIsThrown(String content)
        throws IOException {
      // Arrange
      Path file = write(content);

      // Act & Assert
      assertThrows(RuntimeException.class, () -> reader.read(file));
    }

    @Test
    @DisplayName("given a topology mismatch, when reading, then the file is named in the message")
    void givenTopologyMismatch_whenReading_thenMessageNamesTheFile() throws IOException {
      // Arrange
      Path file = write(">a\nACD\nOO\n>b\nAC\nOO\n");

      // Act
      InvalidInputFileException exception =
          assertThrows(InvalidInputFileException.class, () -> reader.read(file));

      // Assert
      assertEquals(true, exception.getMessage().contains(file.toString()));
    }
  }

  private Path write(String content) throws IOException {
    Path file = directory.resolve("topology.3line");
    Files.writeString(file, content);
    return file;
  }
}
