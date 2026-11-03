package org.tmm2structalign.io;

/**
 * One FASTA entry.
 *
 * @param header header line without the leading {@code >}
 * @param sequence residues, possibly with gaps, without line breaks
 */
public record FastaRecord(String header, String sequence) {

  /** Validates that both fields are present. */
  public FastaRecord {
    if (header == null || sequence == null) {
      throw new IllegalArgumentException("FASTA header and sequence are required");
    }
  }
}
