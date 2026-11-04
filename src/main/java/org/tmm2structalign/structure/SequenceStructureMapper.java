package org.tmm2structalign.structure;

import java.util.Arrays;
import java.util.Locale;
import java.util.Optional;
import org.tmm2structalign.model.ProteinSequence;

/**
 * Maps every residue of an input sequence to a row of its structure's distance matrix, checking
 * that the residues actually agree (review finding A1).
 *
 * <p>Structural row k carries PDB residue number {@code n_k} and, when the PDB sequence is
 * known, residue letter {@code s_k}. Mapping strategies, in the order tried by
 * {@link MappingPolicy#VERIFIED}:
 *
 * <ol>
 *   <li><b>Numbered</b>: residue {@code n_k - 1} of the input maps to row k. Accepted if all
 *       numbers are in range and at least {@link #NUMBERED_MIN_IDENTITY} of the compared
 *       residues are identical (point variants are tolerated).
 *   <li><b>Aligned</b>: semi-global alignment of the PDB and input sequences (end gaps free);
 *       aligned pairs are mapped. Accepted if their identity is at least
 *       {@link #ALIGNED_MIN_IDENTITY}.
 * </ol>
 */
public final class SequenceStructureMapper {

  /** Minimum identity for accepting the PDB numbering. */
  public static final double NUMBERED_MIN_IDENTITY = 0.95;

  /** Minimum identity of aligned pairs for accepting an alignment-based mapping. */
  public static final double ALIGNED_MIN_IDENTITY = 0.90;

  private static final int MATCH = 1;
  private static final int MISMATCH = -1;
  private static final int GAP = -1;

  private static final byte FROM_DIAGONAL = 0;
  private static final byte FROM_UP = 1;
  private static final byte FROM_LEFT = 2;

  private final MappingPolicy policy;

  /**
   * Creates a mapper.
   *
   * @param policy mapping policy
   */
  public SequenceStructureMapper(MappingPolicy policy) {
    this.policy = policy;
  }

  /**
   * Maps a protein onto its structure.
   *
   * @param protein input sequence
   * @param residueNumbers PDB residue number of every structural row
   * @param structureSequence residue letters of the structural rows, if known; when empty the
   *     numbering cannot be checked and the result is {@link MappingStatus#UNVERIFIED}
   * @return the mapping, possibly unusable
   */
  public StructureMapping map(
      ProteinSequence protein, int[] residueNumbers, Optional<String> structureSequence) {
    if (structureSequence.isPresent()
        && structureSequence.get().length() != residueNumbers.length) {
      throw new IllegalArgumentException(
          "Structure sequence of '" + protein.id() + "' has " + structureSequence.get().length()
              + " residues but " + residueNumbers.length + " residue numbers");
    }
    StructureMapping mapping;
    if (structureSequence.isEmpty()) {
      mapping = mapByNumbersUnchecked(protein, residueNumbers);
    } else {
      mapping = mapVerified(protein, residueNumbers, structureSequence.get());
    }
    return mapping;
  }

  private StructureMapping mapVerified(
      ProteinSequence protein, int[] residueNumbers, String structureSequence) {
    StructureMapping numbered = mapByNumbers(protein, residueNumbers, structureSequence);
    StructureMapping mapping = numbered;
    if (!numbered.status().isUsable() && policy == MappingPolicy.VERIFIED) {
      mapping = mapByAlignment(protein, structureSequence, numbered.detail());
    }
    return mapping;
  }

  private static StructureMapping mapByNumbersUnchecked(
      ProteinSequence protein, int[] residueNumbers) {
    StructureMapping mapping;
    Optional<int[]> residueToStructure = numberedMap(protein, residueNumbers);
    if (residueToStructure.isEmpty()) {
      mapping =
          StructureMapping.unusable(
              protein.id(),
              MappingStatus.INVALID_FILES,
              "Residue numbers exceed the sequence length " + protein.length());
    } else {
      int[] map = residueToStructure.get();
      mapping =
          new StructureMapping(
              protein.id(), MappingStatus.UNVERIFIED, 0.0, coverage(map), map,
              "Residue identity not checked");
    }
    return mapping;
  }

  private static StructureMapping mapByNumbers(
      ProteinSequence protein, int[] residueNumbers, String structureSequence) {
    Optional<int[]> residueToStructure = numberedMap(protein, residueNumbers);
    StructureMapping mapping;
    if (residueToStructure.isEmpty()) {
      mapping =
          StructureMapping.unusable(
              protein.id(),
              MappingStatus.LOW_IDENTITY,
              "PDB numbering exceeds the sequence length " + protein.length());
    } else {
      int identical = 0;
      for (int row = 0; row < residueNumbers.length; row++) {
        if (protein.residueAt(residueNumbers[row] - 1) == structureSequence.charAt(row)) {
          identical++;
        }
      }
      double identity = identical / (double) residueNumbers.length;
      int[] map = residueToStructure.get();
      mapping =
          identity >= NUMBERED_MIN_IDENTITY
              ? new StructureMapping(
                  protein.id(), MappingStatus.NUMBERED, identity, coverage(map), map, "")
              : StructureMapping.unusable(
                  protein.id(),
                  MappingStatus.LOW_IDENTITY,
                  String.format(
                      Locale.ROOT, "PDB numbering identity %.3f", identity));
    }
    return mapping;
  }

  private static Optional<int[]> numberedMap(ProteinSequence protein, int[] residueNumbers) {
    int[] residueToStructure = new int[protein.length()];
    Arrays.fill(residueToStructure, -1);
    boolean inRange = true;
    for (int row = 0; row < residueNumbers.length && inRange; row++) {
      int position = residueNumbers[row] - 1;
      inRange = position >= 0 && position < protein.length();
      if (inRange && residueToStructure[position] < 0) {
        residueToStructure[position] = row;
      }
    }
    return inRange ? Optional.of(residueToStructure) : Optional.empty();
  }

  private static StructureMapping mapByAlignment(
      ProteinSequence protein, String structureSequence, String numberingDetail) {
    int[] residueToStructure = semiGlobalAlignment(protein.residues(), structureSequence);
    int aligned = 0;
    int identical = 0;
    for (int residue = 0; residue < residueToStructure.length; residue++) {
      int row = residueToStructure[residue];
      if (row >= 0) {
        aligned++;
        if (protein.residueAt(residue) == structureSequence.charAt(row)) {
          identical++;
        }
      }
    }
    double identity = aligned == 0 ? 0.0 : identical / (double) aligned;
    String detail =
        numberingDetail
            + String.format(Locale.ROOT, "; alignment identity %.3f", identity);
    return identity >= ALIGNED_MIN_IDENTITY
        ? new StructureMapping(
            protein.id(),
            MappingStatus.ALIGNED,
            identity,
            coverage(residueToStructure),
            residueToStructure,
            detail)
        : StructureMapping.unusable(protein.id(), MappingStatus.LOW_IDENTITY, detail);
  }

  /**
   * Aligns two sequences with free end gaps and returns, for every residue of the first, the
   * aligned residue index of the second ({@code -1} if none).
   */
  static int[] semiGlobalAlignment(String sequence, String structure) {
    int rows = sequence.length();
    int columns = structure.length();
    int[] previous = new int[columns + 1];
    int[] current = new int[columns + 1];
    byte[][] trace = new byte[rows + 1][columns + 1];
    for (int column = 1; column <= columns; column++) {
      trace[0][column] = FROM_LEFT;
    }
    int bestScore = Integer.MIN_VALUE;
    int bestRow = rows;
    int bestColumn = columns;
    for (int row = 1; row <= rows; row++) {
      current[0] = 0;
      trace[row][0] = FROM_UP;
      for (int column = 1; column <= columns; column++) {
        int diagonal =
            previous[column - 1]
                + (sequence.charAt(row - 1) == structure.charAt(column - 1) ? MATCH : MISMATCH);
        int up = previous[column] + GAP;
        int left = current[column - 1] + GAP;
        int best = diagonal;
        byte origin = FROM_DIAGONAL;
        if (up > best) {
          best = up;
          origin = FROM_UP;
        }
        if (left > best) {
          best = left;
          origin = FROM_LEFT;
        }
        current[column] = best;
        trace[row][column] = origin;
      }
      if (current[columns] > bestScore) {
        bestScore = current[columns];
        bestRow = row;
        bestColumn = columns;
      }
      int[] swap = previous;
      previous = current;
      current = swap;
    }
    for (int column = 1; column <= columns; column++) {
      if (previous[column] > bestScore) {
        bestScore = previous[column];
        bestRow = rows;
        bestColumn = column;
      }
    }
    return traceBack(trace, rows, bestRow, bestColumn);
  }

  private static int[] traceBack(byte[][] trace, int rows, int endRow, int endColumn) {
    int[] mapping = new int[rows];
    Arrays.fill(mapping, -1);
    int row = endRow;
    int column = endColumn;
    while (row > 0 && column > 0) {
      byte origin = trace[row][column];
      if (origin == FROM_DIAGONAL) {
        mapping[row - 1] = column - 1;
        row--;
        column--;
      } else if (origin == FROM_UP) {
        row--;
      } else {
        column--;
      }
    }
    return mapping;
  }

  private static double coverage(int[] residueToStructure) {
    long mapped = Arrays.stream(residueToStructure).filter(row -> row >= 0).count();
    return residueToStructure.length == 0 ? 0.0 : mapped / (double) residueToStructure.length;
  }
}
