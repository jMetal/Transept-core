package org.transept.model;

/**
 * An unaligned protein sequence with its per-residue predicted topology.
 *
 * @param id identifier, as normalised by the reader that created the sequence
 * @param residues one-letter amino-acid codes, without gaps
 * @param topology one topology code per residue (see {@link Topology#fromCode(char)})
 */
public record ProteinSequence(String id, String residues, String topology) {

  /** Validates that the identifier is present and that every residue has a topology code. */
  public ProteinSequence {
    if (id == null || id.isBlank()) {
      throw new IllegalArgumentException("Protein identifier cannot be blank");
    }
    if (residues == null || residues.isEmpty()) {
      throw new IllegalArgumentException("Protein '" + id + "' has no residues");
    }
    if (residues.indexOf(Alphabet.GAP) >= 0) {
      throw new IllegalArgumentException("Protein '" + id + "' contains gap characters");
    }
    if (topology == null || topology.length() != residues.length()) {
      throw new IllegalArgumentException(
          "Protein '" + id + "' has " + residues.length() + " residues but "
              + (topology == null ? 0 : topology.length()) + " topology codes");
    }
  }

  /**
   * Creates a sequence whose topology is unknown for every residue.
   *
   * @param id identifier
   * @param residues one-letter amino-acid codes, without gaps
   * @return a sequence with {@link Topology#UNKNOWN} topology
   */
  public static ProteinSequence withoutTopology(String id, String residues) {
    String unknown = String.valueOf(Topology.UNKNOWN.code()).repeat(residues.length());
    return new ProteinSequence(id, residues, unknown);
  }

  /** Returns the number of residues. */
  public int length() {
    return residues.length();
  }

  /** Returns the one-letter code of the residue at a 0-based position. */
  public char residueAt(int position) {
    return residues.charAt(position);
  }

  /** Returns the predicted topology of the residue at a 0-based position. */
  public Topology topologyAt(int position) {
    return Topology.fromCode(topology.charAt(position));
  }

  /** Returns true if the residue at a 0-based position lies in a transmembrane helix. */
  public boolean isTransmembrane(int position) {
    return topology.charAt(position) == Topology.TRANSMEMBRANE.code();
  }
}
