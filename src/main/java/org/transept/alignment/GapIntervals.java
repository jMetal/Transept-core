package org.transept.alignment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.transept.model.Alphabet;

/**
 * Gap layout of one aligned row: the aligned length plus the gap runs it contains.
 *
 * <p>Instances are immutable and always canonical: runs are sorted, inclusive, inside
 * {@code [0, alignedLength)}, and separated by at least one residue (adjacent runs are merged).
 * Two rows with the same gaps are therefore {@link #equals(Object) equal}.
 */
public final class GapIntervals {

  private static final int[] EMPTY = new int[0];

  private final int alignedLength;
  private final int[] starts;
  private final int[] ends;

  private GapIntervals(int alignedLength, int[] starts, int[] ends) {
    this.alignedLength = alignedLength;
    this.starts = starts;
    this.ends = ends;
  }

  /**
   * Creates a row without gaps.
   *
   * @param alignedLength number of columns, equal to the number of residues
   * @return a gap-free row
   */
  public static GapIntervals none(int alignedLength) {
    requireNonNegativeLength(alignedLength);
    return new GapIntervals(alignedLength, EMPTY, EMPTY);
  }

  /**
   * Creates a row from inclusive {@code start, end} pairs.
   *
   * @param alignedLength number of columns of the row
   * @param startEndPairs flattened inclusive runs, e.g. {@code 0, 2, 7, 7}
   * @return the canonical row
   * @throws IllegalArgumentException if pairs are incomplete, unsorted, overlapping, adjacent or
   *     out of range
   */
  public static GapIntervals of(int alignedLength, int... startEndPairs) {
    requireNonNegativeLength(alignedLength);
    if (startEndPairs.length % 2 != 0) {
      throw new IllegalArgumentException("Gap runs must be given as start/end pairs");
    }
    int count = startEndPairs.length / 2;
    int[] starts = new int[count];
    int[] ends = new int[count];
    for (int run = 0; run < count; run++) {
      starts[run] = startEndPairs[2 * run];
      ends[run] = startEndPairs[2 * run + 1];
      validateRun(alignedLength, starts, ends, run);
    }
    return new GapIntervals(alignedLength, starts, ends);
  }

  /**
   * Creates a row from a flat list of inclusive {@code start, end} pairs.
   *
   * @param alignedLength number of columns of the row
   * @param startEndPairs flattened inclusive runs
   * @return the canonical row
   */
  public static GapIntervals fromList(int alignedLength, List<Integer> startEndPairs) {
    return of(alignedLength, startEndPairs.stream().mapToInt(Integer::intValue).toArray());
  }

  /**
   * Creates a row from a per-column gap mask, merging consecutive gap columns.
   *
   * @param gapMask {@code true} where the row has a gap
   * @return the canonical row
   */
  public static GapIntervals fromMask(boolean[] gapMask) {
    List<Integer> pairs = new ArrayList<>();
    int column = 0;
    while (column < gapMask.length) {
      if (gapMask[column]) {
        int start = column;
        while (column + 1 < gapMask.length && gapMask[column + 1]) {
          column++;
        }
        pairs.add(start);
        pairs.add(column);
      }
      column++;
    }
    return fromList(gapMask.length, pairs);
  }

  /**
   * Encodes the gaps of an aligned row such as {@code "AC--D-"}.
   *
   * @param alignedRow aligned residues, using {@link Alphabet#GAP} for gaps
   * @return the canonical row
   */
  public static GapIntervals fromAlignedRow(CharSequence alignedRow) {
    boolean[] gapMask = new boolean[alignedRow.length()];
    for (int column = 0; column < gapMask.length; column++) {
      gapMask[column] = alignedRow.charAt(column) == Alphabet.GAP;
    }
    return fromMask(gapMask);
  }

  /** Returns the number of columns of the row. */
  public int alignedLength() {
    return alignedLength;
  }

  /** Returns the number of gap runs. */
  public int runCount() {
    return starts.length;
  }

  /** Returns the first column (inclusive) of a gap run. */
  public int start(int run) {
    return starts[run];
  }

  /** Returns the last column (inclusive) of a gap run. */
  public int end(int run) {
    return ends[run];
  }

  /** Returns the length of a gap run. */
  public int runLength(int run) {
    return ends[run] - starts[run] + 1;
  }

  /** Returns the total number of gap columns. */
  public int gapCount() {
    int gaps = 0;
    for (int run = 0; run < starts.length; run++) {
      gaps += runLength(run);
    }
    return gaps;
  }

  /** Returns the number of residue columns, i.e. the length of the unaligned sequence. */
  public int residueCount() {
    return alignedLength - gapCount();
  }

  /** Returns true if the given column is a gap in this row. */
  public boolean isGapAt(int column) {
    int run = Arrays.binarySearch(starts, column);
    int candidate = run >= 0 ? run : -run - 2;
    return candidate >= 0 && column <= ends[candidate];
  }

  /** Returns a per-column gap mask. */
  public boolean[] toMask() {
    boolean[] mask = new boolean[alignedLength];
    for (int run = 0; run < starts.length; run++) {
      Arrays.fill(mask, starts[run], ends[run] + 1, true);
    }
    return mask;
  }

  /**
   * Maps every column to the 0-based index of the residue it holds.
   *
   * @return array of length {@link #alignedLength()} with {@code -1} at gap columns
   */
  public int[] columnToResidue() {
    int[] residueAtColumn = new int[alignedLength];
    int residue = 0;
    int run = 0;
    for (int column = 0; column < alignedLength; column++) {
      if (run < starts.length && column >= starts[run] && column <= ends[run]) {
        residueAtColumn[column] = -1;
        if (column == ends[run]) {
          run++;
        }
      } else {
        residueAtColumn[column] = residue++;
      }
    }
    return residueAtColumn;
  }

  /**
   * Inserts gaps into an unaligned sequence.
   *
   * @param residues unaligned residues, of length {@link #residueCount()}
   * @return the aligned row
   */
  public String decode(String residues) {
    if (residues.length() != residueCount()) {
      throw new IllegalArgumentException(
          "Sequence has " + residues.length() + " residues but the row expects "
              + residueCount());
    }
    int[] residueAtColumn = columnToResidue();
    StringBuilder aligned = new StringBuilder(alignedLength);
    for (int index : residueAtColumn) {
      aligned.append(index < 0 ? Alphabet.GAP : residues.charAt(index));
    }
    return aligned.toString();
  }

  /** Returns the runs as a flat list of inclusive {@code start, end} pairs. */
  public List<Integer> toList() {
    List<Integer> pairs = new ArrayList<>(2 * starts.length);
    for (int run = 0; run < starts.length; run++) {
      pairs.add(starts[run]);
      pairs.add(ends[run]);
    }
    return List.copyOf(pairs);
  }

  @Override
  public boolean equals(Object other) {
    return other instanceof GapIntervals that
        && alignedLength == that.alignedLength
        && Arrays.equals(starts, that.starts)
        && Arrays.equals(ends, that.ends);
  }

  @Override
  public int hashCode() {
    return 31 * (31 * alignedLength + Arrays.hashCode(starts)) + Arrays.hashCode(ends);
  }

  @Override
  public String toString() {
    return "GapIntervals[length=" + alignedLength + ", runs=" + toList() + "]";
  }

  private static void requireNonNegativeLength(int alignedLength) {
    if (alignedLength < 0) {
      throw new IllegalArgumentException("Aligned length cannot be negative: " + alignedLength);
    }
  }

  private static void validateRun(int alignedLength, int[] starts, int[] ends, int run) {
    if (starts[run] < 0 || ends[run] >= alignedLength || starts[run] > ends[run]) {
      throw new IllegalArgumentException(
          "Invalid gap run [" + starts[run] + ", " + ends[run] + "] for aligned length "
              + alignedLength);
    }
    if (run > 0 && starts[run] <= ends[run - 1] + 1) {
      throw new IllegalArgumentException(
          "Gap runs must be sorted and separated by residues: [" + starts[run - 1] + ", "
              + ends[run - 1] + "] and [" + starts[run] + ", " + ends[run] + "]");
    }
  }
}
