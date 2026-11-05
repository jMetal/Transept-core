package org.transept.testing;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Locale;
import java.util.Map;

/** Writes small input files (PDB, distance matrices) for tests. */
public final class TestFiles {

  private static final Map<Character, String> ONE_TO_THREE =
      Map.ofEntries(
          Map.entry('A', "ALA"),
          Map.entry('R', "ARG"),
          Map.entry('N', "ASN"),
          Map.entry('D', "ASP"),
          Map.entry('C', "CYS"),
          Map.entry('Q', "GLN"),
          Map.entry('E', "GLU"),
          Map.entry('G', "GLY"),
          Map.entry('H', "HIS"),
          Map.entry('I', "ILE"),
          Map.entry('L', "LEU"),
          Map.entry('K', "LYS"),
          Map.entry('M', "MET"),
          Map.entry('F', "PHE"),
          Map.entry('P', "PRO"),
          Map.entry('S', "SER"),
          Map.entry('T', "THR"),
          Map.entry('W', "TRP"),
          Map.entry('Y', "TYR"),
          Map.entry('V', "VAL"));

  private TestFiles() {}

  /** Formats one PDB {@code ATOM} record. */
  public static String atom(String atomName, String residue, char chain, int number) {
    return String.format(
        Locale.ROOT,
        "ATOM  %5d %-4s %3s %c%4d    %8.3f%8.3f%8.3f  1.00 90.00           C%n",
        number,
        " " + atomName,
        residue,
        chain,
        number,
        0.0,
        0.0,
        0.0);
  }

  /** Writes a single-chain PDB file with one Cα per residue, numbered from 1. */
  public static Path writePdb(Path file, String sequence) throws IOException {
    StringBuilder content = new StringBuilder();
    for (int index = 0; index < sequence.length(); index++) {
      content.append(atom("CA", ONE_TO_THREE.get(sequence.charAt(index)), 'A', index + 1));
    }
    Files.writeString(file, content.toString());
    return file;
  }

  /**
   * Writes {@code <id>_D.csv} for residues on a straight line spaced {@code spacing} Å apart, and
   * {@code <id>_idx.csv} numbering them from 1.
   */
  public static void writeLinearStructure(Path directory, String id, int residues, float spacing)
      throws IOException {
    StringBuilder distances = new StringBuilder();
    StringBuilder numbers = new StringBuilder();
    for (int row = 0; row < residues; row++) {
      for (int column = 0; column < residues; column++) {
        distances.append(column == 0 ? "" : ",").append(Math.abs(row - column) * spacing);
      }
      distances.append('\n');
      numbers.append(row == 0 ? "" : ",").append(row + 1);
    }
    Files.writeString(directory.resolve(id + "_D.csv"), distances.toString());
    Files.writeString(directory.resolve(id + "_idx.csv"), numbers.toString());
  }
}
