package org.transept.structure;

import java.util.Optional;

/**
 * Result of mapping one protein onto its structure.
 *
 * @param proteinId protein identifier
 * @param status outcome of the mapping
 * @param identity fraction of identical residues among mapped pairs (0 when not computed)
 * @param coverage fraction of the protein residues with a structural counterpart
 * @param residueToStructure structural index of every residue ({@code -1} if unmapped); empty
 *     array when the structure is not usable
 * @param detail human-readable explanation, empty when there is nothing to report
 */
public record StructureMapping(
    String proteinId,
    MappingStatus status,
    double identity,
    double coverage,
    int[] residueToStructure,
    String detail) {

  /** Creates the result for a protein whose structure is not used. */
  static StructureMapping unusable(String proteinId, MappingStatus status, String detail) {
    return new StructureMapping(proteinId, status, 0.0, 0.0, new int[0], detail);
  }

  /** Returns the explanation, if any. */
  public Optional<String> explanation() {
    return detail.isEmpty() ? Optional.empty() : Optional.of(detail);
  }
}
