package org.transept.alignment;

import java.util.ArrayList;
import java.util.List;

/**
 * Removes columns in which every row has a gap.
 *
 * <p>Replaces the recursive {@code removeGapColumns} of the original implementation, which left
 * some gap-only columns behind (review finding A2).
 */
public final class GapColumnCompactor {

  private GapColumnCompactor() {}

  /**
   * Returns the rows without the columns that are gaps in all of them.
   *
   * @param rows rows of equal aligned length
   * @return compacted rows; the input is returned unchanged when there is nothing to remove
   * @throws IllegalArgumentException if rows are empty or have different lengths
   */
  public static List<GapIntervals> compact(List<GapIntervals> rows) {
    boolean[] occupied = occupiedColumns(rows);
    int[] newIndex = new int[occupied.length];
    int kept = 0;
    for (int column = 0; column < occupied.length; column++) {
      newIndex[column] = kept;
      if (occupied[column]) {
        kept++;
      }
    }

    List<GapIntervals> compacted = rows;
    if (kept != occupied.length) {
      compacted = new ArrayList<>(rows.size());
      for (GapIntervals row : rows) {
        compacted.add(remap(row, occupied, newIndex, kept));
      }
      compacted = List.copyOf(compacted);
    }
    return compacted;
  }

  /**
   * Counts the columns in which every row has a gap.
   *
   * @param rows rows of equal aligned length
   * @return number of gap-only columns
   */
  public static int countGapOnlyColumns(List<GapIntervals> rows) {
    int gapOnly = 0;
    for (boolean columnOccupied : occupiedColumns(rows)) {
      if (!columnOccupied) {
        gapOnly++;
      }
    }
    return gapOnly;
  }

  private static boolean[] occupiedColumns(List<GapIntervals> rows) {
    if (rows.isEmpty()) {
      throw new IllegalArgumentException("Cannot compact an empty alignment");
    }
    int length = rows.get(0).alignedLength();
    boolean[] occupied = new boolean[length];
    for (GapIntervals row : rows) {
      if (row.alignedLength() != length) {
        throw new IllegalArgumentException(
            "Rows have different lengths: " + length + " and " + row.alignedLength());
      }
      boolean[] gaps = row.toMask();
      for (int column = 0; column < length; column++) {
        occupied[column] |= !gaps[column];
      }
    }
    return occupied;
  }

  private static GapIntervals remap(
      GapIntervals row, boolean[] occupied, int[] newIndex, int newLength) {
    boolean[] gaps = row.toMask();
    boolean[] newGaps = new boolean[newLength];
    for (int column = 0; column < gaps.length; column++) {
      if (occupied[column] && gaps[column]) {
        newGaps[newIndex[column]] = true;
      }
    }
    return GapIntervals.fromMask(newGaps);
  }
}
