package org.transept.structure;

/** Outcome of mapping a protein sequence onto its structure. */
public enum MappingStatus {
  /** PDB residue numbers index the sequence and the residues agree. */
  NUMBERED(true),
  /** Numbering disagreed; the mapping comes from aligning the PDB and input sequences. */
  ALIGNED(true),
  /** Residue numbers used without checking residues (no PDB file to check them). */
  UNVERIFIED(true),
  /** Neither numbering nor alignment reached the identity threshold; structure not used. */
  LOW_IDENTITY(false),
  /** Distance or residue-number files are missing; structure not used. */
  MISSING_FILES(false),
  /** Structural files exist but are malformed; structure not used. */
  INVALID_FILES(false);

  private final boolean usable;

  MappingStatus(boolean usable) {
    this.usable = usable;
  }

  /** Returns true if the structure is used by structural objectives. */
  public boolean isUsable() {
    return usable;
  }
}
