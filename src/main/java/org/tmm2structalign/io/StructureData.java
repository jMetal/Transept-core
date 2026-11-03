package org.tmm2structalign.io;

/**
 * Cα distance matrix of a structure and the residue number of each matrix row.
 *
 * <p>Arrays are shared, not copied, because distance matrices can be large; treat them as
 * read-only.
 *
 * @param distances symmetric matrix of Cα distances in Å, one row per structural residue
 * @param residueNumbers residue number (PDB {@code resSeq}) of each row
 */
public record StructureData(float[][] distances, int[] residueNumbers) {

  /** Validates that every matrix row has a residue number. */
  public StructureData {
    if (distances == null || residueNumbers == null || distances.length != residueNumbers.length) {
      throw new IllegalArgumentException("Distance matrix and residue numbers must match");
    }
  }

  /** Returns the number of structural residues. */
  public int residueCount() {
    return residueNumbers.length;
  }
}
