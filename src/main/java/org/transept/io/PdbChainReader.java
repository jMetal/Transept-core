package org.transept.io;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Extracts the residue sequence of one chain from a PDB file.
 *
 * <p>Follows the conventions of {@code distanceandindex.py}, which built the distance matrices:
 * first model only, chain {@code A} if present (otherwise the first chain), standard residues
 * ({@code ATOM} records) that have a Cα atom, first alternate location. The k-th residue
 * returned therefore corresponds to row k of the distance matrix.
 */
public final class PdbChainReader {

  private static final String PREFERRED_CHAIN = "A";

  private static final Map<String, Character> THREE_TO_ONE =
      Map.ofEntries(
          Map.entry("ALA", 'A'),
          Map.entry("ARG", 'R'),
          Map.entry("ASN", 'N'),
          Map.entry("ASP", 'D'),
          Map.entry("CYS", 'C'),
          Map.entry("GLN", 'Q'),
          Map.entry("GLU", 'E'),
          Map.entry("GLY", 'G'),
          Map.entry("HIS", 'H'),
          Map.entry("ILE", 'I'),
          Map.entry("LEU", 'L'),
          Map.entry("LYS", 'K'),
          Map.entry("MET", 'M'),
          Map.entry("PHE", 'F'),
          Map.entry("PRO", 'P'),
          Map.entry("SER", 'S'),
          Map.entry("THR", 'T'),
          Map.entry("TRP", 'W'),
          Map.entry("TYR", 'Y'),
          Map.entry("VAL", 'V'),
          Map.entry("SEC", 'U'),
          Map.entry("PYL", 'O'));

  /**
   * Residue sequence of a chain and the residue number of each residue.
   *
   * @param sequence one-letter residue codes ({@code X} for non-standard residues)
   * @param residueNumbers PDB residue numbers, one per residue
   */
  public record PdbChain(String sequence, int[] residueNumbers) {}

  /**
   * Reads the chain sequence of a PDB file.
   *
   * @param file PDB file
   * @return the chain sequence
   * @throws IOException if the file cannot be read
   * @throws InvalidInputFileException if the file has no Cα atoms
   */
  public PdbChain read(Path file) throws IOException {
    List<String> caLines = new ArrayList<>();
    for (String line : Files.readAllLines(file, StandardCharsets.UTF_8)) {
      if (line.startsWith("ENDMDL")) {
        break;
      }
      if (line.startsWith("ATOM") && line.length() >= 27 && "CA".equals(atomName(line))) {
        caLines.add(line);
      }
    }
    if (caLines.isEmpty()) {
      throw new InvalidInputFileException("No C-alpha atoms in PDB file " + file);
    }
    String chain = selectChain(caLines);
    StringBuilder sequence = new StringBuilder();
    List<Integer> numbers = new ArrayList<>();
    Set<String> seenResidues = new HashSet<>();
    for (String line : caLines) {
      String residueKey = line.substring(22, 27);
      if (chainOf(line).equals(chain) && seenResidues.add(residueKey)) {
        sequence.append(THREE_TO_ONE.getOrDefault(line.substring(17, 20).trim(), 'X'));
        numbers.add(Integer.parseInt(line.substring(22, 26).trim()));
      }
    }
    return new PdbChain(
        sequence.toString(), numbers.stream().mapToInt(Integer::intValue).toArray());
  }

  private static String atomName(String line) {
    return line.substring(12, 16).trim();
  }

  private static String chainOf(String line) {
    return line.substring(21, 22);
  }

  private static String selectChain(List<String> caLines) {
    return caLines.stream()
        .map(PdbChainReader::chainOf)
        .filter(PREFERRED_CHAIN::equals)
        .findFirst()
        .orElse(chainOf(caLines.get(0)));
  }
}
