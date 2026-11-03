package org.tmm2structalign.testing;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import org.tmm2structalign.alignment.GapIntervals;
import org.tmm2structalign.alignment.MsaSolution;
import org.tmm2structalign.model.ProteinSequence;
import org.tmm2structalign.model.ProteinSet;
import org.tmm2structalign.structure.MappingStatus;
import org.tmm2structalign.structure.ProteinStructures;
import org.tmm2structalign.structure.StructureMapping;

/** Builders of small proteins, solutions and structures for tests. */
public final class Alignments {

  private Alignments() {}

  /** Creates proteins named p0, p1, ... with the given residues and topologies. */
  public static ProteinSet proteins(String... residuesAndTopologies) {
    List<ProteinSequence> proteins = new ArrayList<>();
    for (int index = 0; index < residuesAndTopologies.length; index++) {
      String[] parts = residuesAndTopologies[index].split(":");
      proteins.add(new ProteinSequence("p" + index, parts[0], parts[1]));
    }
    return new ProteinSet(proteins);
  }

  /** Creates a two-objective solution from aligned rows. */
  public static MsaSolution solution(ProteinSet proteins, String... alignedRows) {
    List<GapIntervals> rows = new ArrayList<>();
    for (String row : alignedRows) {
      rows.add(GapIntervals.fromAlignedRow(row));
    }
    return new MsaSolution(proteins, rows, 2);
  }

  /** Returns the aligned rows of a solution. */
  public static List<String> rows(MsaSolution solution) {
    List<String> rows = new ArrayList<>();
    for (int sequence = 0; sequence < solution.variables().size(); sequence++) {
      rows.add(solution.toMsa().alignedRow(sequence));
    }
    return rows;
  }

  /** Creates a random valid solution with up to {@code maxGaps} gaps inserted per row. */
  public static MsaSolution randomSolution(ProteinSet proteins, int maxGaps, Random random) {
    int longest = proteins.proteins().stream().mapToInt(ProteinSequence::length).max().orElse(0);
    int length = longest + maxGaps;
    List<GapIntervals> rows = new ArrayList<>();
    for (ProteinSequence protein : proteins.proteins()) {
      boolean[] gaps = new boolean[length];
      int toPlace = length - protein.length();
      while (toPlace > 0) {
        int column = random.nextInt(length);
        if (!gaps[column]) {
          gaps[column] = true;
          toPlace--;
        }
      }
      rows.add(GapIntervals.fromMask(gaps));
    }
    return new MsaSolution(proteins, rows, 2);
  }

  /** Structures in which every residue of every protein is mapped to its own row. */
  public static ProteinStructures fullyStructured(ProteinSet proteins) {
    List<StructureMapping> mappings = new ArrayList<>();
    List<float[][]> distances = new ArrayList<>();
    for (ProteinSequence protein : proteins.proteins()) {
      int[] identity = new int[protein.length()];
      for (int residue = 0; residue < identity.length; residue++) {
        identity[residue] = residue;
      }
      mappings.add(
          new StructureMapping(protein.id(), MappingStatus.NUMBERED, 1, 1, identity, ""));
      distances.add(new float[protein.length()][protein.length()]);
    }
    return new ProteinStructures(mappings, distances);
  }

  /** Structures in which no protein has structural data. */
  public static ProteinStructures unstructured(ProteinSet proteins) {
    List<StructureMapping> mappings = new ArrayList<>();
    List<float[][]> distances = new ArrayList<>();
    for (ProteinSequence protein : proteins.proteins()) {
      mappings.add(
          new StructureMapping(protein.id(), MappingStatus.MISSING_FILES, 0, 0, new int[0], ""));
      distances.add(null);
    }
    return new ProteinStructures(mappings, distances);
  }

  /**
   * Counts the gap runs of a row placed between two consecutive transmembrane residues, i.e.
   * the runs that break a helix.
   */
  public static int helixBreakingRuns(ProteinSequence protein, GapIntervals row) {
    int[] residueAtColumn = row.columnToResidue();
    int breaking = 0;
    int residuesBefore = 0;
    int run = 0;
    for (int column = 0; column < residueAtColumn.length && run < row.runCount(); column++) {
      if (column == row.start(run)) {
        if (residuesBefore > 0
            && residuesBefore < protein.length()
            && protein.isTransmembrane(residuesBefore - 1)
            && protein.isTransmembrane(residuesBefore)) {
          breaking++;
        }
        run++;
      }
      residuesBefore += residueAtColumn[column] >= 0 ? 1 : 0;
    }
    return breaking;
  }
}
