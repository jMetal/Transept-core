package org.tmm2structalign.example;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.function.Supplier;
import org.tmm2structalign.alignment.GapIntervals;
import org.tmm2structalign.alignment.Msa;
import org.tmm2structalign.alignment.MsaSolution;
import org.tmm2structalign.model.ProteinSequence;
import org.tmm2structalign.model.ProteinSet;
import org.tmm2structalign.objective.LddtAggregation;
import org.tmm2structalign.objective.Objective;
import org.tmm2structalign.objective.StructuralLddt;
import org.tmm2structalign.objective.TopologyAwareSumOfPairs;
import org.tmm2structalign.operator.MergeAdjacentGapRunsMutation;
import org.tmm2structalign.operator.RandomSource;
import org.tmm2structalign.operator.TmSafeGapShiftMutation;
import org.tmm2structalign.operator.TmSafeSinglePointCrossover;
import org.tmm2structalign.problem.StructuralMsaProblem;
import org.tmm2structalign.structure.MappingStatus;
import org.tmm2structalign.structure.ProteinStructures;
import org.tmm2structalign.structure.StructureMapping;
import org.uma.jmetal.operator.mutation.MutationOperator;
import org.uma.jmetal.util.pseudorandom.JMetalRandom;

/**
 * Companion of the tutorial "T2. Problem formulation, encoding and variation operators".
 *
 * <p>Builds three short proteins with one transmembrane helix each, synthetic structures, and a
 * seed alignment, and prints the encoding, the objectives and the effect of every variation
 * operator. Run it with
 *
 * <pre>
 * java -cp target/TM-MSAStructAligner-1.0-jar-with-dependencies.jar \
 *     org.tmm2structalign.example.FormulationTour
 * </pre>
 *
 * <p>The output is deterministic, and a test checks that the tutorial shows it unchanged.
 */
public final class FormulationTour {

  /** Seed of the global jMetal generator used by the operators. */
  public static final long RANDOM_SEED = 1L;

  /** Children shown for every mutation. */
  private static final int CHILDREN = 4;

  private static final String[][] PROTEINS = {
    {"p1", "PQTLLVAVLGNKE", "OOOMMMMMMMIII"},
    {"p2", "PTLLIAVLGKE", "OOMMMMMMMII"},
    {"p3", "PQSTLVAILGNE", "OOOOMMMMMMII"},
  };

  /** Alignment from which the synthetic structures are built: its columns are homologous. */
  private static final List<String> SEED =
      List.of("PQTLLVAVLGNKE", "P-TLLIAVLG-KE", "PQSTLVAILGN-E");

  private final ProteinSet proteins;
  private final ProteinStructures structures;
  private final StructuralMsaProblem problem;

  private FormulationTour() {
    List<ProteinSequence> sequences = new ArrayList<>();
    for (String[] protein : PROTEINS) {
      sequences.add(new ProteinSequence(protein[0], protein[1], protein[2]));
    }
    proteins = new ProteinSet(sequences);
    structures = syntheticStructures(proteins);
    List<Objective> objectives =
        List.of(
            TopologyAwareSumOfPairs.withDefaultParameters(),
            new StructuralLddt(
                structures, StructuralLddt.DEFAULT_CUTOFF, LddtAggregation.HARMONIC_OVER_LENGTH));
    problem = new StructuralMsaProblem("tour", proteins, objectives);
  }

  /**
   * Prints the tour.
   *
   * @param args not used
   */
  public static void main(String[] args) {
    System.out.print(run());
  }

  /**
   * Runs the tour.
   *
   * @return the text printed by {@link #main(String[])}
   */
  public static String run() {
    return new FormulationTour().tour();
  }

  private String tour() {
    StringBuilder out = new StringBuilder();
    proteinsSection(out);
    encodingSection(out);
    objectivesSection(out);
    compactionSection(out);
    crossoverSection(out);
    mutationSection(out);
    return out.toString();
  }

  private void proteinsSection(StringBuilder out) {
    heading(out, "1. Proteins");
    for (ProteinSequence protein : proteins.proteins()) {
      out.append(String.format(Locale.ROOT, "%s  %s%n    %s%n", protein.id(),
          protein.residues(), protein.topology()));
    }
  }

  private void encodingSection(StringBuilder out) {
    heading(out, "2. Alignment and encoding");
    MsaSolution seed = solution(SEED);
    for (int sequence = 0; sequence < proteins.size(); sequence++) {
      GapIntervals row = seed.variables().get(sequence);
      out.append(String.format(Locale.ROOT, "%s  %s  gap runs %s%n",
          proteins.get(sequence).id(), SEED.get(sequence), runs(row)));
    }
    out.append(String.format(Locale.ROOT, "length %d%n", seed.alignmentLength()));
    out.append("transmembrane columns (M), gaps included:\n");
    Msa msa = seed.toMsa();
    for (int sequence = 0; sequence < msa.size(); sequence++) {
      StringBuilder flags = new StringBuilder();
      for (int column = 0; column < msa.length(); column++) {
        flags.append(msa.isTransmembrane(sequence, column) ? 'M' : '.');
      }
      out.append(String.format(Locale.ROOT, "%s  %s%n", proteins.get(sequence).id(), flags));
    }
  }

  private void objectivesSection(StringBuilder out) {
    heading(out, "3. Objectives");
    describe(out, "seed alignment", solution(SEED));
    describe(out, "p2 shifted one column",
        solution(List.of("PQTLLVAVLGNKE-", "P--TLLIAVLGKE-", "PQSTLVAILGN-E-")));
  }

  private void compactionSection(StringBuilder out) {
    heading(out, "4. Gap-only columns");
    MsaSolution padded =
        solution(List.of("PQTLLVAVLG-NKE", "P-TLLIAVLG--KE", "PQSTLVAILG-N-E"));
    out.append("before evaluation:\n");
    rows(out, padded);
    problem.evaluate(padded);
    out.append("after evaluation:\n");
    rows(out, padded);
  }

  private void crossoverSection(StringBuilder out) {
    heading(out, "5. Crossover");
    MsaSolution first = solution(SEED);
    out.append("first parent:\n");
    rows(out, first);
    MsaSolution second = solution(List.of("PQTLLVAVLGNKE--", "--PTLLIAVLGKE--", "PQSTLVAILGNE---"));
    out.append("second parent:\n");
    rows(out, second);
    TmSafeSinglePointCrossover single = new TmSafeSinglePointCrossover(1.0);
    out.append(String.format(Locale.ROOT, "valid cuts of the first parent: %s%n",
        single.compatibleCuts(first)));
    out.append("single-point child, cut at column 1:\n");
    describe(out, null, single.recombine(first, second, 1));
  }

  private void mutationSection(StringBuilder out) {
    heading(out, "6. Mutations");
    MsaSolution parent = solution(List.of("PQT-LLVAVLGNKE", "P-TLLIAVLG-KE-", "PQSTLVAILGN--E"));
    out.append("parent:\n");
    describe(out, null, parent);
    RandomSource random = RandomSource.jmetal();
    Map<String, Supplier<MutationOperator<MsaSolution>>> mutations = new LinkedHashMap<>();
    mutations.put("TM_SAFE_GAP_SHIFT", () ->
        new TmSafeGapShiftMutation(1.0, structures, 60, random));
    mutations.put("MERGE_GAP_RUNS", () -> new MergeAdjacentGapRunsMutation(1.0, random));
    for (Map.Entry<String, Supplier<MutationOperator<MsaSolution>>> entry : mutations.entrySet()) {
      out.append(String.format(Locale.ROOT, "%s, %d children:%n", entry.getKey(), CHILDREN));
      JMetalRandom.getInstance().setSeed(RANDOM_SEED);
      MutationOperator<MsaSolution> mutation = entry.getValue().get();
      for (int child = 0; child < CHILDREN; child++) {
        describe(out, null, mutation.execute(parent.copy()));
      }
    }
  }

  private void describe(StringBuilder out, String title, MsaSolution solution) {
    problem.evaluate(solution);
    if (title != null) {
      out.append(title).append(":\n");
    }
    rows(out, solution);
    out.append(String.format(Locale.ROOT, "    SOPwTP %.1f   lDDT %.4f   broken helices %d%n",
        problem.naturalValue(solution, 0), problem.naturalValue(solution, 1),
        brokenHelices(solution)));
  }

  /** Counts the gap runs placed between two consecutive transmembrane residues. */
  private static int brokenHelices(MsaSolution solution) {
    int broken = 0;
    for (int sequence = 0; sequence < solution.variables().size(); sequence++) {
      ProteinSequence protein = solution.proteins().get(sequence);
      GapIntervals row = solution.variables().get(sequence);
      int[] residueAtColumn = row.columnToResidue();
      for (int run = 0; run < row.runCount(); run++) {
        int before = 0;
        for (int column = 0; column < row.start(run); column++) {
          before += residueAtColumn[column] >= 0 ? 1 : 0;
        }
        if (before > 0 && before < protein.length() && protein.isTransmembrane(before - 1)
            && protein.isTransmembrane(before)) {
          broken++;
        }
      }
    }
    return broken;
  }

  private static void rows(StringBuilder out, MsaSolution solution) {
    for (int sequence = 0; sequence < solution.variables().size(); sequence++) {
      out.append("    ").append(solution.toMsa().alignedRow(sequence)).append('\n');
    }
  }

  private static String runs(GapIntervals row) {
    List<String> runs = new ArrayList<>();
    for (int run = 0; run < row.runCount(); run++) {
      runs.add(row.start(run) == row.end(run)
          ? Integer.toString(row.start(run))
          : row.start(run) + "-" + row.end(run));
    }
    return runs.toString();
  }

  private static void heading(StringBuilder out, String title) {
    out.append("== ").append(title).append(" ==\n");
  }

  private MsaSolution solution(List<String> rows) {
    return problem.createSolution(Msa.fromAlignedRows(proteins, rows));
  }

  /**
   * Builds structures in which residues aligned in the same column of {@link #SEED} sit at
   * nearly the same position of an ideal helix, so the seed alignment is structurally good.
   */
  private static ProteinStructures syntheticStructures(ProteinSet proteins) {
    List<StructureMapping> mappings = new ArrayList<>();
    List<float[][]> distances = new ArrayList<>();
    for (int sequence = 0; sequence < proteins.size(); sequence++) {
      ProteinSequence protein = proteins.get(sequence);
      double[][] coordinates = new double[protein.length()][];
      int residue = 0;
      String row = SEED.get(sequence);
      for (int column = 0; column < row.length(); column++) {
        if (row.charAt(column) != '-') {
          double angle = Math.toRadians(100.0 * column);
          double stretch = 1.0 + 0.03 * sequence;
          coordinates[residue++] =
              new double[] {2.3 * Math.cos(angle), 2.3 * Math.sin(angle), 1.5 * column * stretch};
        }
      }
      float[][] matrix = new float[protein.length()][protein.length()];
      int[] identity = new int[protein.length()];
      for (int i = 0; i < protein.length(); i++) {
        identity[i] = i;
        for (int j = 0; j < protein.length(); j++) {
          double dx = coordinates[i][0] - coordinates[j][0];
          double dy = coordinates[i][1] - coordinates[j][1];
          double dz = coordinates[i][2] - coordinates[j][2];
          matrix[i][j] = (float) Math.sqrt(dx * dx + dy * dy + dz * dz);
        }
      }
      mappings.add(
          new StructureMapping(protein.id(), MappingStatus.NUMBERED, 1.0, 1.0, identity, ""));
      distances.add(matrix);
    }
    return new ProteinStructures(mappings, distances);
  }
}
