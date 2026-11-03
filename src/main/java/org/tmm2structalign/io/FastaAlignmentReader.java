package org.tmm2structalign.io;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import org.tmm2structalign.alignment.Msa;
import org.tmm2structalign.model.Alphabet;
import org.tmm2structalign.model.ProteinSet;

/**
 * Reads aligned FASTA files and matches their entries to a known set of proteins.
 *
 * <p>Entries may appear in any order and with tool-specific identifiers (see
 * {@link ProteinIds}). The alignment is returned in the canonical order of the protein set.
 * Residues are upper-cased and {@code .} is read as a gap.
 */
public final class FastaAlignmentReader {

  /**
   * Reads every entry of a FASTA file.
   *
   * @param file path of the file
   * @return entries in file order
   * @throws IOException if the file cannot be read
   * @throws InvalidInputFileException if sequence lines appear before the first header
   */
  public List<FastaRecord> readRecords(Path file) throws IOException {
    List<FastaRecord> records = new ArrayList<>();
    String header = null;
    StringBuilder sequence = new StringBuilder();
    for (String rawLine : Files.readAllLines(file, StandardCharsets.UTF_8)) {
      String line = rawLine.trim();
      if (line.startsWith(">")) {
        if (header != null) {
          records.add(new FastaRecord(header, sequence.toString()));
        }
        header = line.substring(1).trim();
        sequence.setLength(0);
      } else if (!line.isEmpty()) {
        if (header == null) {
          throw new InvalidInputFileException("Sequence data before the first header in " + file);
        }
        sequence.append(line);
      }
    }
    if (header != null) {
      records.add(new FastaRecord(header, sequence.toString()));
    }
    return records;
  }

  /**
   * Reads an aligned FASTA file and returns it in the canonical order of the proteins.
   *
   * @param file path of the aligned FASTA file
   * @param proteins expected proteins; every one must appear exactly once
   * @return the alignment
   * @throws IOException if the file cannot be read
   * @throws InvalidInputFileException if proteins are missing, duplicated, unexpected, have
   *     different aligned lengths or do not match their sequences
   */
  public Msa read(Path file, ProteinSet proteins) throws IOException {
    Map<String, String> alignedById = new HashMap<>();
    for (FastaRecord record : readRecords(file)) {
      String id = ProteinIds.normalize(record.header());
      String aligned = record.sequence().toUpperCase(Locale.ROOT).replace('.', Alphabet.GAP);
      if (alignedById.put(id, aligned) != null) {
        throw new InvalidInputFileException("Duplicated protein '" + id + "' in " + file);
      }
    }

    List<String> rows = new ArrayList<>(proteins.size());
    for (int index = 0; index < proteins.size(); index++) {
      String id = ProteinIds.normalize(proteins.get(index).id());
      String aligned = alignedById.remove(id);
      if (aligned == null) {
        throw new InvalidInputFileException("Protein '" + id + "' not found in " + file);
      }
      if (!rows.isEmpty() && aligned.length() != rows.get(0).length()) {
        throw new InvalidInputFileException(
            "Aligned length of '" + id + "' is " + aligned.length() + " but expected "
                + rows.get(0).length() + " in " + file);
      }
      rows.add(aligned);
    }
    if (!alignedById.isEmpty()) {
      throw new InvalidInputFileException(
          "Unexpected proteins " + alignedById.keySet() + " in " + file);
    }

    try {
      return Msa.fromAlignedRows(proteins, rows);
    } catch (IllegalArgumentException exception) {
      throw new InvalidInputFileException(
          "Invalid alignment in " + file + ": " + exception.getMessage(), exception);
    }
  }
}
