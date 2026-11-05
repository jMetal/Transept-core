package org.transept.io;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.ValueSource;
import org.transept.alignment.Msa;
import org.transept.model.ProteinSequence;
import org.transept.model.ProteinSet;

@DisplayName("Unit tests for FASTA input and output: ProteinIds, FastaAlignmentReader/Writer")
class FastaIoTest {

  @TempDir Path directory;

  private FastaAlignmentReader reader;
  private ProteinSet proteins;

  @BeforeEach
  void setUp() {
    reader = new FastaAlignmentReader();
    proteins =
        new ProteinSet(
            List.of(
                ProteinSequence.withoutTopology("ada2b_human", "ACDE"),
                ProteinSequence.withoutTopology("P12345", "FGH")));
  }

  @Nested
  @DisplayName("When normalising identifiers")
  class Identifiers {

    @ParameterizedTest
    @CsvSource({
      "'>ada2b_human', ada2b_human",
      "'ADA2B_HUMAN description text', ada2b_human",
      "'ada2b_human.pdb', ada2b_human",
      "'ada2b_human.pdb:A', ada2b_human",
      "'ada2b_human:B', ada2b_human",
      "'sp|P12345|ADA2B_HUMAN', p12345",
      "'tr|Q99999|X', q99999",
      "'ada2b_human|SP+TM', ada2b_human"
    })
    @DisplayName("given a tool-specific header, when normalising, then the id is canonical")
    void givenToolSpecificHeader_whenNormalising_thenIdIsCanonical(String header, String id) {
      // Act & Assert
      assertEquals(id, ProteinIds.normalize(header));
    }

    @ParameterizedTest
    @ValueSource(strings = {">", "  ", ".pdb", ":A"})
    @DisplayName("given a header without identifier, when normalising, then an exception is thrown")
    void givenEmptyIdentifier_whenNormalising_thenInvalidInputFileExceptionIsThrown(
        String header) {
      // Act & Assert
      assertThrows(InvalidInputFileException.class, () -> ProteinIds.normalize(header));
    }

    @Test
    @DisplayName("given a null header, when normalising, then an exception is thrown")
    void givenNullHeader_whenNormalising_thenInvalidInputFileExceptionIsThrown() {
      // Act & Assert
      assertThrows(InvalidInputFileException.class, () -> ProteinIds.normalize(null));
    }
  }

  @Nested
  @DisplayName("When reading aligned FASTA files")
  class Reading {

    @Test
    @DisplayName("given shuffled tool output, when reading, then the canonical order is restored")
    void givenShuffledToolOutput_whenReading_thenCanonicalOrderIsRestored() throws IOException {
      // Arrange
      Path file = write(">sp|P12345|X\nf.g-\nH-\n>ada2b_human.pdb:A\nAC\nDE--\n");

      // Act
      Msa msa = reader.read(file, proteins);

      // Assert
      assertEquals("ACDE--", msa.alignedRow(0));
      assertEquals("F-G-H-", msa.alignedRow(1));
    }

    @ParameterizedTest
    @ValueSource(
        strings = {
          ">ada2b_human\nACDE\n", // missing protein
          ">ada2b_human\nACDE\n>p12345\nFGH-\n>extra\nAAAA\n", // unexpected protein
          ">ada2b_human\nACDE\n>ada2b_human.pdb\nACDE\n>p12345\nFGH-\n", // duplicate
          ">ada2b_human\nACDE\n>p12345\nFGH--\n", // different lengths
          ">ada2b_human\nACDW\n>p12345\nFGH-\n", // different residues
          "ACDE\n>ada2b_human\nACDE\n" // data before the first header
        })
    @DisplayName("given an inconsistent file, when reading, then an exception is thrown")
    void givenInconsistentFile_whenReading_thenInvalidInputFileExceptionIsThrown(String content)
        throws IOException {
      // Arrange
      Path file = write(content);

      // Act & Assert
      assertThrows(InvalidInputFileException.class, () -> reader.read(file, proteins));
    }
  }

  @Nested
  @DisplayName("When writing aligned FASTA files")
  class Writing {

    @Test
    @DisplayName("given an alignment, when writing and reading back, then it is unchanged")
    void givenAlignment_whenWritingAndReadingBack_thenItIsUnchanged() throws IOException {
      // Arrange
      Msa msa = Msa.fromAlignedRows(proteins, List.of("AC-DE", "-FGH-"));
      Path file = directory.resolve("out.fasta");

      // Act
      new FastaAlignmentWriter(2).write(msa, file);
      Msa readBack = reader.read(file, proteins);

      // Assert
      assertEquals(msa.rows(), readBack.rows());
      assertTrue(Files.readString(file).startsWith(">ada2b_human\nAC\n-D\nE\n>P12345\n"));
    }

    @Test
    @DisplayName("given a non-positive line width, when creating a writer, then it is rejected")
    void givenNonPositiveLineWidth_whenCreatingWriter_thenIllegalArgumentExceptionIsThrown() {
      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> new FastaAlignmentWriter(0));
    }

    @Test
    @DisplayName("given the default writer, when formatting, then lines hold 60 residues")
    void givenDefaultWriter_whenFormatting_thenLinesHoldSixtyResidues() {
      // Arrange
      String longResidues = "A".repeat(61);
      ProteinSet longProteins =
          new ProteinSet(
              List.of(
                  ProteinSequence.withoutTopology("a", longResidues),
                  ProteinSequence.withoutTopology("b", longResidues)));
      Msa msa = Msa.fromAlignedRows(longProteins, List.of(longResidues, longResidues));

      // Act
      String text = new FastaAlignmentWriter().format(msa);

      // Assert
      assertEquals(60, text.split("\n")[1].length());
      assertEquals(1, text.split("\n")[2].length());
    }
  }

  private Path write(String content) throws IOException {
    Path file = directory.resolve("aligned.fasta");
    Files.writeString(file, content);
    return file;
  }
}
