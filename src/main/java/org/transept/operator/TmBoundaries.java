package org.transept.operator;

import org.transept.model.ProteinSequence;

/**
 * Helpers shared by the operators that never break a transmembrane helix.
 *
 * <p>Rows are handled as per-column gap masks. An <b>insertion position</b> {@code k} of a row of
 * length {@code L} is a value in {@code [0, L]}: inserting at {@code k} places the new gaps
 * before column {@code k}, after the residues that precede it. The position is valid when those
 * residues and the next one are not both transmembrane.
 */
final class TmBoundaries {

  private TmBoundaries() {}

  /**
   * Tells whether gaps placed after the given number of residues separate two transmembrane
   * residues.
   *
   * @param protein the protein of the row
   * @param residues residues before the boundary
   * @return true if the boundary lies inside a helix
   */
  static boolean splitsHelix(ProteinSequence protein, int residues) {
    return residues > 0
        && residues < protein.length()
        && protein.isTransmembrane(residues - 1)
        && protein.isTransmembrane(residues);
  }

  /**
   * Counts the residues before every insertion position of a row.
   *
   * @param gaps gap mask of the row
   * @return array of length {@code gaps.length + 1}
   */
  static int[] residuesBefore(boolean[] gaps) {
    int[] counts = new int[gaps.length + 1];
    for (int column = 0; column < gaps.length; column++) {
      counts[column + 1] = counts[column] + (gaps[column] ? 0 : 1);
    }
    return counts;
  }

  /**
   * Removes consecutive columns from a mask.
   *
   * @param mask gap mask
   * @param start first column to remove
   * @param count number of columns to remove
   * @return a new, shorter mask
   */
  static boolean[] removeColumns(boolean[] mask, int start, int count) {
    boolean[] result = new boolean[mask.length - count];
    System.arraycopy(mask, 0, result, 0, start);
    System.arraycopy(mask, start + count, result, start, mask.length - start - count);
    return result;
  }

  /**
   * Inserts gap columns into a mask.
   *
   * @param mask gap mask
   * @param position insertion position, in {@code [0, mask.length]}
   * @param count number of gaps to insert
   * @return a new, longer mask
   */
  static boolean[] insertGaps(boolean[] mask, int position, int count) {
    boolean[] result = new boolean[mask.length + count];
    System.arraycopy(mask, 0, result, 0, position);
    for (int index = position; index < position + count; index++) {
      result[index] = true;
    }
    System.arraycopy(mask, position, result, position + count, mask.length - position);
    return result;
  }
}
