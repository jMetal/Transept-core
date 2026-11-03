package org.tmm2structalign.alignment;

import java.util.ArrayList;
import java.util.List;
import org.tmm2structalign.model.Alphabet;
import org.tmm2structalign.model.ProteinSequence;
import org.tmm2structalign.model.ProteinSet;

/**
 * Immutable multiple sequence alignment decoded into integer lookup tables.
 *
 * <p>Objectives read alignments through this class so that a solution is decoded once per
 * evaluation, without allocating per-residue objects.
 *
 * <p>Transmembrane status of a gap follows the original TM-M2Align convention: a gap counts as
 * transmembrane only when both residues flanking the gap run are transmembrane. At the ends of a
 * sequence the terminal residue is used for both flanks.
 */
public final class Msa {

  private final ProteinSet proteins;
  private final List<GapIntervals> rows;
  private final int length;
  private final int[][] residueAtColumn;
  private final boolean[][] transmembrane;

  private Msa(ProteinSet proteins, List<GapIntervals> rows) {
    this.proteins = proteins;
    this.rows = List.copyOf(rows);
    this.length = rows.get(0).alignedLength();
    this.residueAtColumn = new int[rows.size()][];
    this.transmembrane = new boolean[rows.size()][];
    for (int sequence = 0; sequence < rows.size(); sequence++) {
      residueAtColumn[sequence] = rows.get(sequence).columnToResidue();
      transmembrane[sequence] =
          transmembraneFlags(proteins.get(sequence), residueAtColumn[sequence]);
    }
  }

  /**
   * Creates an alignment from one gap layout per protein.
   *
   * @param proteins proteins in canonical order
   * @param rows gap layout of each protein, in the same order
   * @return the alignment
   * @throws IllegalArgumentException if the rows do not match the proteins or have different
   *     lengths
   */
  public static Msa of(ProteinSet proteins, List<GapIntervals> rows) {
    if (rows.size() != proteins.size()) {
      throw new IllegalArgumentException(
          "Expected " + proteins.size() + " rows but got " + rows.size());
    }
    int length = rows.get(0).alignedLength();
    for (int sequence = 0; sequence < rows.size(); sequence++) {
      GapIntervals row = rows.get(sequence);
      ProteinSequence protein = proteins.get(sequence);
      if (row.alignedLength() != length) {
        throw new IllegalArgumentException(
            "Row of '" + protein.id() + "' has length " + row.alignedLength()
                + " but the alignment has length " + length);
      }
      if (row.residueCount() != protein.length()) {
        throw new IllegalArgumentException(
            "Row of '" + protein.id() + "' holds " + row.residueCount()
                + " residues but the protein has " + protein.length());
      }
    }
    return new Msa(proteins, rows);
  }

  /**
   * Creates an alignment from aligned strings, checking that each one spells its protein.
   *
   * @param proteins proteins in canonical order
   * @param alignedRows aligned strings in the same order; residue case is ignored
   * @return the alignment
   * @throws IllegalArgumentException if a row does not contain exactly its protein's residues
   */
  public static Msa fromAlignedRows(ProteinSet proteins, List<String> alignedRows) {
    if (alignedRows.size() != proteins.size()) {
      throw new IllegalArgumentException(
          "Expected " + proteins.size() + " aligned rows but got " + alignedRows.size());
    }
    List<GapIntervals> rows = new ArrayList<>(alignedRows.size());
    for (int sequence = 0; sequence < alignedRows.size(); sequence++) {
      String aligned = alignedRows.get(sequence);
      ProteinSequence protein = proteins.get(sequence);
      String ungapped = aligned.replace(String.valueOf(Alphabet.GAP), "");
      if (!ungapped.equalsIgnoreCase(protein.residues())) {
        throw new IllegalArgumentException(
            "Aligned row of '" + protein.id() + "' does not match its sequence (expected "
                + protein.length() + " residues, found " + ungapped.length() + ")");
      }
      rows.add(GapIntervals.fromAlignedRow(aligned));
    }
    return of(proteins, rows);
  }

  /** Returns the proteins in canonical order. */
  public ProteinSet proteins() {
    return proteins;
  }

  /** Returns the gap layout of every row. */
  public List<GapIntervals> rows() {
    return rows;
  }

  /** Returns the number of sequences. */
  public int size() {
    return rows.size();
  }

  /** Returns the number of columns. */
  public int length() {
    return length;
  }

  /** Returns the 0-based residue index at a column, or {@code -1} for a gap. */
  public int residueIndex(int sequence, int column) {
    return residueAtColumn[sequence][column];
  }

  /** Returns true if the sequence has a gap at the column. */
  public boolean isGap(int sequence, int column) {
    return residueAtColumn[sequence][column] < 0;
  }

  /** Returns the residue letter at a column, or {@link Alphabet#GAP}. */
  public char letter(int sequence, int column) {
    int residue = residueAtColumn[sequence][column];
    return residue < 0 ? Alphabet.GAP : proteins.get(sequence).residueAt(residue);
  }

  /** Returns the transmembrane status of a column, including gaps (see class comment). */
  public boolean isTransmembrane(int sequence, int column) {
    return transmembrane[sequence][column];
  }

  /** Returns the number of sequences with a residue at the column. */
  public int occupancy(int column) {
    int residues = 0;
    for (int[] row : residueAtColumn) {
      if (row[column] >= 0) {
        residues++;
      }
    }
    return residues;
  }

  /** Returns the aligned string of a sequence. */
  public String alignedRow(int sequence) {
    return rows.get(sequence).decode(proteins.get(sequence).residues());
  }

  private static boolean[] transmembraneFlags(ProteinSequence protein, int[] residueAtColumn) {
    boolean[] flags = new boolean[residueAtColumn.length];
    int lastResidue = protein.length() - 1;
    int residuesSeen = 0;
    for (int column = 0; column < residueAtColumn.length; column++) {
      int residue = residueAtColumn[column];
      if (residue >= 0) {
        flags[column] = protein.isTransmembrane(residue);
        residuesSeen++;
      } else {
        int left = Math.max(residuesSeen - 1, 0);
        int right = Math.min(residuesSeen, lastResidue);
        flags[column] = protein.isTransmembrane(left) && protein.isTransmembrane(right);
      }
    }
    return flags;
  }
}
