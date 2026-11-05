package org.transept.testing;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * Writes a small dataset with the directory layout of {@code Dataset_Structural_GPCRs/GPCRdb}:
 * five related proteins with two transmembrane segments, helix-like Cα coordinates, PDB files
 * and three seed alignments ({@code left}, {@code right} and {@code true}).
 */
public final class SyntheticDataset {

  /** Dataset name used in file names. */
  public static final String NAME = "synthetic";

  /** Seed methods available in {@code precomputed/synthetic/}. */
  public static final List<String> SEED_METHODS = List.of("left", "right", "true");

  private static final String ANCESTOR =
      "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ";
  private static final String TOPOLOGY =
      "OOOOOMMMMMMMMMMMMMMMMMMMMIIIIIIIIIIMMMMMMMMMMMMMMMMMMMMOOOOOOOOOOO";
  private static final int[][] DELETIONS = {{}, {3}, {60, 61}, {1, 2, 64}, {40}};

  private SyntheticDataset() {}

  /**
   * Writes the dataset under a root directory.
   *
   * @param root directory that plays the role of {@code GPCRdb}
   * @return the root
   * @throws IOException if files cannot be written
   */
  public static Path write(Path root) throws IOException {
    List<String> trueRows = new ArrayList<>();
    StringBuilder topology = new StringBuilder();
    Path distances = Files.createDirectories(root.resolve("distances").resolve(NAME));
    Path pdbs = Files.createDirectories(root.resolve("alphafold_pdb").resolve(NAME));
    List<String> sequences = new ArrayList<>();
    for (int protein = 0; protein < DELETIONS.length; protein++) {
      StringBuilder residues = new StringBuilder();
      StringBuilder labels = new StringBuilder();
      StringBuilder aligned = new StringBuilder();
      for (int position = 0; position < ANCESTOR.length(); position++) {
        if (contains(DELETIONS[protein], position)) {
          aligned.append('-');
        } else {
          residues.append(ANCESTOR.charAt(position));
          labels.append(TOPOLOGY.charAt(position));
          aligned.append(ANCESTOR.charAt(position));
        }
      }
      String id = "prot" + protein + "_synth";
      topology.append('>').append(id).append(" | TM\n").append(residues).append('\n')
          .append(labels).append('\n');
      trueRows.add(aligned.toString());
      sequences.add(residues.toString());
      writeStructure(distances, pdbs, id, residues.toString(), DELETIONS[protein]);
    }
    Path topologyDirectory =
        Files.createDirectories(root.resolve("sequences").resolve("tmregions"));
    Files.writeString(
        topologyDirectory.resolve(NAME + "_predicted_topologies.3line"), topology.toString());

    Path precomputed = Files.createDirectories(root.resolve("precomputed").resolve(NAME));
    writeFasta(precomputed.resolve("true_" + NAME + ".fasta"), trueRows);
    writeFasta(precomputed.resolve("left_" + NAME + ".fasta"), pad(sequences, true));
    writeFasta(precomputed.resolve("right_" + NAME + ".fasta"), pad(sequences, false));
    return root;
  }

  private static void writeStructure(
      Path distances, Path pdbs, String id, String residues, int[] deletions) throws IOException {
    List<double[]> coordinates = new ArrayList<>();
    for (int position = 0; position < ANCESTOR.length(); position++) {
      if (!contains(deletions, position)) {
        double angle = Math.toRadians(100.0 * position);
        coordinates.add(
            new double[] {2.3 * Math.cos(angle), 2.3 * Math.sin(angle), 1.5 * position});
      }
    }
    StringBuilder matrix = new StringBuilder();
    StringBuilder numbers = new StringBuilder();
    for (int row = 0; row < coordinates.size(); row++) {
      for (int column = 0; column < coordinates.size(); column++) {
        matrix.append(column == 0 ? "" : ",")
            .append(String.format(Locale.ROOT, "%.3f", distance(coordinates, row, column)));
      }
      matrix.append('\n');
      numbers.append(row == 0 ? "" : ",").append(row + 1);
    }
    Files.writeString(distances.resolve(id + "_D.csv"), matrix.toString());
    Files.writeString(distances.resolve(id + "_idx.csv"), numbers.toString());
    TestFiles.writePdb(pdbs.resolve(id + ".pdb"), residues);
  }

  private static double distance(List<double[]> coordinates, int first, int second) {
    double[] a = coordinates.get(first);
    double[] b = coordinates.get(second);
    double dx = a[0] - b[0];
    double dy = a[1] - b[1];
    double dz = a[2] - b[2];
    return Math.sqrt(dx * dx + dy * dy + dz * dz);
  }

  private static List<String> pad(List<String> sequences, boolean left) {
    int longest = sequences.stream().mapToInt(String::length).max().orElse(0);
    return sequences.stream()
        .map(sequence -> {
          String gaps = "-".repeat(longest - sequence.length());
          return left ? sequence + gaps : gaps + sequence;
        })
        .toList();
  }

  private static void writeFasta(Path file, List<String> rows) throws IOException {
    StringBuilder fasta = new StringBuilder();
    for (int protein = 0; protein < rows.size(); protein++) {
      fasta.append(">prot").append(protein).append("_synth\n").append(rows.get(protein))
          .append('\n');
    }
    Files.writeString(file, fasta.toString());
  }

  private static boolean contains(int[] values, int value) {
    for (int candidate : values) {
      if (candidate == value) {
        return true;
      }
    }
    return false;
  }
}
