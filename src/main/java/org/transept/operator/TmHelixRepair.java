package org.transept.operator;

import java.util.ArrayList;
import java.util.List;
import org.transept.alignment.GapIntervals;
import org.transept.alignment.Msa;
import org.transept.model.ProteinSequence;

/**
 * Removes the gaps that split a transmembrane helix from an alignment, for instance a seed
 * alignment produced by a tool that does not know the topology.
 *
 * <p>In every row, each gap run placed between two consecutive transmembrane residues is moved to
 * the nearest boundary that does not split a helix, measured in residues (ties go towards the
 * N-terminus), and this is repeated until no run splits a helix. Moving a run of {@code k} gaps
 * across {@code d} residues shifts those residues {@code k} columns; the length of the row and of
 * the alignment does not change. The residues of every row are kept in order.
 */
public final class TmHelixRepair {

  /**
   * A repaired alignment.
   *
   * @param alignment the alignment without gaps inside helices
   * @param movedRuns gap runs that were moved
   */
  public record Result(Msa alignment, int movedRuns) {}

  private TmHelixRepair() {}

  /**
   * Repairs every row of an alignment.
   *
   * @param msa alignment to repair
   * @return the repaired alignment and the number of gap runs moved
   */
  public static Result repair(Msa msa) {
    List<GapIntervals> rows = new ArrayList<>(msa.size());
    int moved = 0;
    for (int sequence = 0; sequence < msa.size(); sequence++) {
      ProteinSequence protein = msa.proteins().get(sequence);
      boolean[] mask = msa.rows().get(sequence).toMask();
      int[] start = new int[1];
      while (findSplittingRun(protein, mask, start)) {
        mask = moveRun(protein, mask, start[0]);
        moved++;
      }
      rows.add(GapIntervals.fromMask(mask));
    }
    return new Result(Msa.of(msa.proteins(), rows), moved);
  }

  /** Finds the first gap run of a row that splits a helix; its first column goes to start[0]. */
  private static boolean findSplittingRun(ProteinSequence protein, boolean[] mask, int[] start) {
    int[] before = TmBoundaries.residuesBefore(mask);
    for (int column = 0; column < mask.length; column++) {
      boolean runStarts = mask[column] && (column == 0 || !mask[column - 1]);
      if (runStarts && TmBoundaries.splitsHelix(protein, before[column])) {
        start[0] = column;
        return true;
      }
    }
    return false;
  }

  /** Moves the gap run that starts at a column to the nearest boundary outside the helices. */
  private static boolean[] moveRun(ProteinSequence protein, boolean[] mask, int start) {
    int end = start;
    while (end + 1 < mask.length && mask[end + 1]) {
      end++;
    }
    int length = end - start + 1;
    int residues = TmBoundaries.residuesBefore(mask)[start];
    int target = nearestSafeBoundary(protein, residues);
    boolean[] removed = TmBoundaries.removeColumns(mask, start, length);
    return TmBoundaries.insertGaps(removed, positionAfter(removed, target), length);
  }

  /** Number of residues before the nearest boundary that does not split a helix. */
  private static int nearestSafeBoundary(ProteinSequence protein, int residues) {
    for (int distance = 1; ; distance++) {
      if (residues - distance >= 0 && !TmBoundaries.splitsHelix(protein, residues - distance)) {
        return residues - distance;
      }
      if (residues + distance <= protein.length()
          && !TmBoundaries.splitsHelix(protein, residues + distance)) {
        return residues + distance;
      }
    }
  }

  /** Insertion position right after the given number of residues (before the first residue). */
  private static int positionAfter(boolean[] mask, int residues) {
    if (residues == 0) {
      int column = 0;
      while (column < mask.length && mask[column]) {
        column++;
      }
      return column;
    }
    int seen = 0;
    for (int column = 0; column < mask.length; column++) {
      if (!mask[column] && ++seen == residues) {
        return column + 1;
      }
    }
    return mask.length;
  }
}
