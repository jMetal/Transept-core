package org.transept.structure;

/** How residue numbers in structural files are related to positions in the input sequence. */
public enum MappingPolicy {
  /**
   * Default. Use the PDB numbering when the residues agree; otherwise align the PDB sequence to
   * the input sequence; reject the structure if both fall below the identity thresholds.
   */
  VERIFIED,
  /** Use the PDB numbering only when the residues agree; never fall back to alignment. */
  NUMBERED_ONLY
}
