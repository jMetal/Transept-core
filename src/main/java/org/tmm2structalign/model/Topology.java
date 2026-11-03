package org.tmm2structalign.model;

/**
 * Predicted membrane topology of a residue, as encoded in the third line of a 3-line topology
 * record.
 */
public enum Topology {
  TRANSMEMBRANE('M'),
  INSIDE('I'),
  OUTSIDE('O'),
  SIGNAL_PEPTIDE('S'),
  UNKNOWN('U');

  private final char code;

  Topology(char code) {
    this.code = code;
  }

  /** Returns the one-letter code used in topology files. */
  public char code() {
    return code;
  }

  /** Returns true for residues inside a transmembrane helix. */
  public boolean isTransmembrane() {
    return this == TRANSMEMBRANE;
  }

  /**
   * Maps a one-letter topology code to its value. Unrecognised codes map to {@link #UNKNOWN}.
   *
   * @param code topology letter, case-sensitive
   * @return the corresponding topology
   */
  public static Topology fromCode(char code) {
    return switch (code) {
      case 'M' -> TRANSMEMBRANE;
      case 'I' -> INSIDE;
      case 'O' -> OUTSIDE;
      case 'S' -> SIGNAL_PEPTIDE;
      default -> UNKNOWN;
    };
  }
}
