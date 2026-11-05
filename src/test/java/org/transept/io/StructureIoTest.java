package org.transept.io;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
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
import org.junit.jupiter.params.provider.CsvSource;
import org.transept.io.PdbChainReader.PdbChain;

@DisplayName("Unit tests for structural input: StructureDataReader and PdbChainReader")
class StructureIoTest {

  @TempDir Path directory;

  @Nested
  @DisplayName("When reading distance and residue-number files")
  class DistanceFiles {

    private StructureDataReader reader;

    @BeforeEach
    void setUp() {
      reader = new StructureDataReader();
    }

    @Test
    @DisplayName("given valid files, when reading, then distances and numbers are returned")
    void givenValidFiles_whenReading_thenDistancesAndNumbersAreReturned() throws IOException {
      // Arrange
      Path distances = write("p_D.csv", "0,3.8,7\n3.8,0,3.8\n7,3.8,0\n");
      Path index = write("p_idx.csv", "1,2\n3");

      // Act
      StructureData data = reader.read(distances, index);

      // Assert
      assertEquals(3, data.residueCount());
      assertEquals(7.0f, data.distances()[0][2]);
      assertArrayEquals(new int[] {1, 2, 3}, data.residueNumbers());
      assertEquals(directory.resolve("x_D.csv"), StructureDataReader.distanceFile(directory, "x"));
      assertEquals(directory.resolve("x_idx.csv"), StructureDataReader.indexFile(directory, "x"));
    }

    @ParameterizedTest
    @CsvSource(
        delimiter = '|',
        value = {
          "0,1\\n1,0,2|1,2", // not square
          "0,1\\n2,0|1,2", // not symmetric
          "1,1\\n1,0|1,2", // non-zero diagonal
          "0,1\\n1,0|1,2,3", // fewer rows than numbers
          "0,1\\n1,0|2,1", // decreasing numbers
          "0,1\\n1,0|0,1", // non-positive number
          "0,x\\nx,0|1,2", // not a number
          "0,1\\n1,0|1,b", // not an integer
          "' '|1" // empty matrix
        })
    @DisplayName("given inconsistent files, when reading, then an exception is thrown")
    void givenInconsistentFiles_whenReading_thenInvalidInputFileExceptionIsThrown(
        String distanceContent, String indexContent) throws IOException {
      // Arrange
      Path distances = write("q_D.csv", distanceContent.replace("\\n", "\n"));
      Path index = write("q_idx.csv", indexContent);

      // Act & Assert
      assertThrows(InvalidInputFileException.class, () -> reader.read(distances, index));
    }

    @Test
    @DisplayName("given an empty residue-number file, when reading, then an exception is thrown")
    void givenEmptyIndexFile_whenReading_thenInvalidInputFileExceptionIsThrown()
        throws IOException {
      // Arrange
      Path distances = write("r_D.csv", "0\n");
      Path index = write("r_idx.csv", "\n");

      // Act & Assert
      assertThrows(InvalidInputFileException.class, () -> reader.read(distances, index));
    }
  }

  @Nested
  @DisplayName("When reading PDB chains")
  class PdbFiles {

    private PdbChainReader reader;

    @BeforeEach
    void setUp() {
      reader = new PdbChainReader();
    }

    @Test
    @DisplayName("given a two-model file, when reading, then chain A of model 1 is returned")
    void givenTwoModels_whenReading_thenPreferredChainOfFirstModelIsReturned() throws IOException {
      // Arrange
      Path pdb =
          write(
              "x.pdb",
              atom("N", "MET", 'A', 1)
                  + atom("CA", "MET", 'A', 1)
                  + atom("CA", "GLY", 'A', 2)
                  + atom("CA", "GLY", 'A', 2) // alternate location, ignored
                  + atom("CA", "MSE", 'A', 3) // non-standard residue
                  + atom("CA", "TRP", 'B', 4) // other chain
                  + "HETATM    1  O   HOH A 100      0.000   0.000   0.000\n"
                  + "ENDMDL\n"
                  + atom("CA", "LYS", 'A', 9));

      // Act
      PdbChain chain = reader.read(pdb);

      // Assert
      assertEquals("MGX", chain.sequence());
      assertArrayEquals(new int[] {1, 2, 3}, chain.residueNumbers());
    }

    @Test
    @DisplayName("given no chain A, when reading, then the first chain is used")
    void givenNoChainA_whenReading_thenFirstChainIsUsed() throws IOException {
      // Arrange
      Path pdb = write("y.pdb", atom("CA", "ALA", 'B', 5) + atom("CA", "CYS", 'C', 6));

      // Act & Assert
      assertEquals("A", reader.read(pdb).sequence());
    }

    @Test
    @DisplayName("given no C-alpha atoms, when reading, then an exception is thrown")
    void givenNoCalphaAtoms_whenReading_thenInvalidInputFileExceptionIsThrown()
        throws IOException {
      // Arrange
      Path pdb = write("z.pdb", atom("N", "ALA", 'A', 1));

      // Act & Assert
      assertThrows(InvalidInputFileException.class, () -> reader.read(pdb));
    }
  }

  private Path write(String name, String content) throws IOException {
    Path file = directory.resolve(name);
    Files.writeString(file, content);
    return file;
  }

  private static String atom(String atomName, String residue, char chain, int number) {
    return String.format(
        "ATOM  %5d %-4s %3s %c%4d    %8.3f%8.3f%8.3f  1.00 90.00           C%n",
        number, " " + atomName, residue, chain, number, 0.0, 0.0, 0.0);
  }
}
