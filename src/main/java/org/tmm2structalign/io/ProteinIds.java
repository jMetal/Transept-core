package org.tmm2structalign.io;

import java.util.Locale;

/**
 * Normalises protein identifiers found in FASTA headers so that the same protein matches across
 * files produced by different tools.
 *
 * <p>Rules, applied in order: drop a leading {@code >}, keep the first whitespace-separated
 * token, take the accession of UniProt headers ({@code sp|P12345|NAME} gives {@code P12345}) or
 * the first {@code |}-separated field otherwise, drop a chain suffix ({@code :A}), drop a
 * {@code .pdb} extension, and lower-case the result.
 */
public final class ProteinIds {

  private ProteinIds() {}

  /**
   * Returns the normalised identifier of a FASTA header.
   *
   * @param header header line, with or without the leading {@code >}
   * @return the normalised, lower-case identifier
   * @throws InvalidInputFileException if the header yields an empty identifier
   */
  public static String normalize(String header) {
    if (header == null) {
      throw new InvalidInputFileException("FASTA header is missing");
    }
    String value = header.trim();
    if (value.startsWith(">")) {
      value = value.substring(1).trim();
    }
    value = value.isEmpty() ? value : value.split("\\s+", 2)[0];

    String[] fields = value.split("\\|", -1);
    boolean uniProt =
        fields.length >= 2
            && ("sp".equalsIgnoreCase(fields[0]) || "tr".equalsIgnoreCase(fields[0]));
    value = uniProt ? fields[1] : fields[0];

    int chainSeparator = value.indexOf(':');
    if (chainSeparator >= 0) {
      value = value.substring(0, chainSeparator);
    }
    value = value.trim().toLowerCase(Locale.ROOT);
    if (value.endsWith(".pdb")) {
      value = value.substring(0, value.length() - ".pdb".length()).trim();
    }
    if (value.isEmpty()) {
      throw new InvalidInputFileException("Empty protein identifier in header: " + header);
    }
    return value;
  }
}
