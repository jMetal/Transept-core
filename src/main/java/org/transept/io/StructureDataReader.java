package org.transept.io;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

/**
 * Reads the {@code <id>_D.csv} distance matrix and the {@code <id>_idx.csv} residue numbers
 * produced by {@code distanceandindex.py} in the dataset repository.
 */
public final class StructureDataReader {

  private static final float SYMMETRY_TOLERANCE = 1.0e-3f;

  /**
   * Returns the distance-matrix path of a protein.
   *
   * @param directory directory holding the structural files of one dataset
   * @param proteinId identifier as written in the topology file
   * @return the {@code _D.csv} path
   */
  public static Path distanceFile(Path directory, String proteinId) {
    return directory.resolve(proteinId + "_D.csv");
  }

  /**
   * Returns the residue-number path of a protein.
   *
   * @param directory directory holding the structural files of one dataset
   * @param proteinId identifier as written in the topology file
   * @return the {@code _idx.csv} path
   */
  public static Path indexFile(Path directory, String proteinId) {
    return directory.resolve(proteinId + "_idx.csv");
  }

  /**
   * Reads and validates the structural files of a protein.
   *
   * @param distanceFile {@code _D.csv} file: square, symmetric, zero diagonal
   * @param indexFile {@code _idx.csv} file: positive, non-decreasing residue numbers
   * @return the structure data
   * @throws IOException if a file cannot be read
   * @throws InvalidInputFileException if the content is malformed or inconsistent
   */
  public StructureData read(Path distanceFile, Path indexFile) throws IOException {
    float[][] distances = readDistances(distanceFile);
    int[] residueNumbers = readResidueNumbers(indexFile);
    if (distances.length != residueNumbers.length) {
      throw new InvalidInputFileException(
          "Distance matrix " + distanceFile + " has " + distances.length + " rows but "
              + indexFile + " has " + residueNumbers.length + " residue numbers");
    }
    validateResidueNumbers(indexFile, residueNumbers);
    validateSymmetry(distanceFile, distances);
    return new StructureData(distances, residueNumbers);
  }

  private static float[][] readDistances(Path file) throws IOException {
    List<String> rows =
        Files.readAllLines(file, StandardCharsets.UTF_8).stream()
            .map(String::trim)
            .filter(line -> !line.isEmpty())
            .toList();
    if (rows.isEmpty()) {
      throw new InvalidInputFileException("Empty distance matrix: " + file);
    }
    int size = rows.size();
    float[][] distances = new float[size][];
    for (int row = 0; row < size; row++) {
      String[] values = rows.get(row).split(",");
      if (values.length != size) {
        throw new InvalidInputFileException(
            "Distance matrix " + file + " is not square: row " + row + " has " + values.length
                + " values, expected " + size);
      }
      distances[row] = new float[size];
      for (int column = 0; column < size; column++) {
        distances[row][column] = parseFloat(file, values[column]);
      }
    }
    return distances;
  }

  private static int[] readResidueNumbers(Path file) throws IOException {
    String content = Files.readString(file, StandardCharsets.UTF_8).trim();
    if (content.isEmpty()) {
      throw new InvalidInputFileException("Empty residue-number file: " + file);
    }
    String[] values = content.split("[,\\s]+");
    int[] numbers = new int[values.length];
    for (int index = 0; index < values.length; index++) {
      try {
        numbers[index] = Integer.parseInt(values[index]);
      } catch (NumberFormatException exception) {
        throw new InvalidInputFileException(
            "Invalid residue number '" + values[index] + "' in " + file, exception);
      }
    }
    return numbers;
  }

  private static void validateResidueNumbers(Path file, int[] numbers) {
    for (int index = 0; index < numbers.length; index++) {
      if (numbers[index] <= 0 || (index > 0 && numbers[index] < numbers[index - 1])) {
        throw new InvalidInputFileException(
            "Residue numbers in " + file + " must be positive and non-decreasing (position "
                + index + ", value " + numbers[index] + ")");
      }
    }
  }

  private static void validateSymmetry(Path file, float[][] distances) {
    for (int row = 0; row < distances.length; row++) {
      if (Math.abs(distances[row][row]) > SYMMETRY_TOLERANCE) {
        throw new InvalidInputFileException(
            "Non-zero diagonal at " + row + " in distance matrix " + file);
      }
      for (int column = row + 1; column < distances.length; column++) {
        if (Math.abs(distances[row][column] - distances[column][row]) > SYMMETRY_TOLERANCE) {
          throw new InvalidInputFileException(
              "Distance matrix " + file + " is not symmetric at (" + row + ", " + column + ")");
        }
      }
    }
  }

  private static float parseFloat(Path file, String value) {
    try {
      return Float.parseFloat(value.trim());
    } catch (NumberFormatException exception) {
      throw new InvalidInputFileException(
          "Invalid distance '" + value + "' in " + file, exception);
    }
  }
}
