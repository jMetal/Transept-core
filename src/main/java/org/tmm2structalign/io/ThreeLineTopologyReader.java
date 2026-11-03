package org.tmm2structalign.io;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import org.tmm2structalign.model.ProteinSequence;
import org.tmm2structalign.model.ProteinSet;

/**
 * Reads 3-line topology files: a header {@code >id | annotation}, the sequence, and one
 * topology letter per residue. Blank lines between records are ignored.
 *
 * <p>The protein identifier is the header text before the first {@code |}, trimmed, as in the
 * original implementation; it is used verbatim to locate structural files.
 */
public final class ThreeLineTopologyReader {

  /**
   * Reads a 3-line topology file.
   *
   * @param file path of the file
   * @return the proteins in file order
   * @throws IOException if the file cannot be read
   * @throws InvalidInputFileException if a record is malformed
   */
  public ProteinSet read(Path file) throws IOException {
    List<String> lines =
        Files.readAllLines(file, StandardCharsets.UTF_8).stream()
            .map(String::trim)
            .filter(line -> !line.isEmpty())
            .toList();
    if (lines.size() % 3 != 0) {
      throw new InvalidInputFileException(
          "Incomplete 3-line record in " + file + ": " + lines.size() + " non-empty lines");
    }
    List<ProteinSequence> proteins = new ArrayList<>(lines.size() / 3);
    for (int first = 0; first < lines.size(); first += 3) {
      proteins.add(parseRecord(file, lines.get(first), lines.get(first + 1), lines.get(first + 2)));
    }
    return new ProteinSet(proteins);
  }

  private static ProteinSequence parseRecord(
      Path file, String header, String residues, String topology) {
    if (!header.startsWith(">")) {
      throw new InvalidInputFileException(
          "Expected a header starting with '>' in " + file + " but found: " + header);
    }
    int separator = header.indexOf('|');
    String id = (separator > 0 ? header.substring(1, separator) : header.substring(1)).trim();
    try {
      return new ProteinSequence(id, residues.toUpperCase(Locale.ROOT), topology);
    } catch (IllegalArgumentException exception) {
      throw new InvalidInputFileException(
          "Invalid record in " + file + ": " + exception.getMessage(), exception);
    }
  }
}
