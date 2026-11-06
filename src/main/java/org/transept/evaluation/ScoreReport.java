package org.transept.evaluation;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.PathMatcher;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;

/**
 * Writes a CSV report that compares every alignment of a directory tree with a reference alignment.
 *
 * <p>An alignment that cannot be compared does not stop the report: its row has the status {@code
 * failed} and explains why.
 */
public final class ScoreReport {

  private static final int METRIC_COLUMNS = 16;

  private ScoreReport() {}

  /**
   * How many alignments a report compared.
   *
   * @param scored alignments compared with the reference
   * @param failed alignments that could not be compared
   */
  public record Summary(int scored, int failed) {}

  /**
   * Compares every alignment of a directory tree with a reference alignment.
   *
   * @param reference reference alignment (FASTA)
   * @param topology three-line topology file of the dataset
   * @param root directory searched recursively for alignments
   * @param glob pattern of the names of the alignment files
   * @param output CSV report to write
   * @return how many alignments were compared
   * @throws IOException if the inputs cannot be read, no alignment matches, or the report cannot be
   *     written
   */
  public static Summary write(Path reference, Path topology, Path root, String glob, Path output)
      throws IOException {
    List<Path> alignments = discoverAlignments(root, glob);
    if (alignments.isEmpty()) {
      throw new IOException("No files matching '" + glob + "' were found under " + root);
    }
    validateOutputTarget(output, reference, topology, alignments);
    AlignmentEvaluator evaluator = new AlignmentEvaluator(reference, topology);
    Path outputParent = output.toAbsolutePath().getParent();
    if (outputParent != null) {
      Files.createDirectories(outputParent);
    }
    int failed = 0;
    try (BufferedWriter writer = Files.newBufferedWriter(output, StandardCharsets.UTF_8)) {
      writer.write(csvHeader());
      writer.newLine();
      for (Path alignment : alignments) {
        String row;
        try {
          row = scoredRow(root, alignment, evaluator.evaluate(alignment));
        } catch (IOException | IllegalArgumentException exception) {
          failed++;
          row = failedRow(root, alignment, exception.getMessage());
        }
        writer.write(row);
        writer.newLine();
      }
    }
    return new Summary(alignments.size() - failed, failed);
  }

  private static void validateOutputTarget(
      Path output, Path reference, Path topology, List<Path> alignments) throws IOException {
    Path target = output.toAbsolutePath().normalize();
    if (target.equals(reference.toAbsolutePath().normalize())
        || target.equals(topology.toAbsolutePath().normalize())
        || alignments.stream().anyMatch(path -> target.equals(path.toAbsolutePath().normalize()))) {
      throw new IOException("The report must not overwrite an input file: " + output);
    }
  }

  static List<Path> discoverAlignments(Path root, String glob) throws IOException {
    PathMatcher matcher = FileSystems.getDefault().getPathMatcher("glob:" + glob);
    List<Path> matches = new ArrayList<>();
    try (var paths = Files.walk(root)) {
      paths
          .filter(Files::isRegularFile)
          .filter(path -> matcher.matches(path.getFileName()))
          .forEach(matches::add);
    }
    matches.sort(Comparator.comparing(path -> root.relativize(path).toString()));
    return matches;
  }

  private static String csvHeader() {
    return String.join(
        ",",
        "relative_path",
        "run_id",
        "msa_file",
        "sequence_count",
        "pair_precision",
        "pair_recall",
        "pair_f1",
        "tm_pair_precision",
        "tm_pair_recall",
        "tm_pair_f1",
        "tm_gap_rate",
        "reference_pairs",
        "test_pairs",
        "correct_pairs",
        "reference_tm_pairs",
        "test_tm_pairs",
        "correct_tm_pairs",
        "tm_gap_count",
        "tm_gap_opportunities",
        "status",
        "message");
  }

  private static String location(Path root, Path alignment) {
    Path relative = root.relativize(alignment);
    String runId =
        relative.getNameCount() > 1 ? relative.getParent().getFileName().toString() : "";
    return String.join(
        ",", csv(relative.toString()), csv(runId), csv(alignment.getFileName().toString()));
  }

  private static String scoredRow(Path root, Path alignment, AlignmentEvaluationResult result) {
    return String.join(
        ",",
        location(root, alignment),
        Integer.toString(result.sequenceCount()),
        decimal(result.pairPrecision()),
        decimal(result.pairRecall()),
        decimal(result.pairF1()),
        decimal(result.tmPairPrecision()),
        decimal(result.tmPairRecall()),
        decimal(result.tmPairF1()),
        decimal(result.tmGapRate()),
        Long.toString(result.referencePairCount()),
        Long.toString(result.testPairCount()),
        Long.toString(result.correctPairCount()),
        Long.toString(result.referenceTmPairCount()),
        Long.toString(result.testTmPairCount()),
        Long.toString(result.correctTmPairCount()),
        Long.toString(result.tmGapCount()),
        Long.toString(result.tmGapOpportunityCount()),
        "ok",
        csv(""));
  }

  private static String failedRow(Path root, Path alignment, String message) {
    List<String> cells = new ArrayList<>(List.of(location(root, alignment)));
    cells.addAll(Collections.nCopies(METRIC_COLUMNS, ""));
    cells.add("failed");
    cells.add(csv(message == null ? "" : message));
    return String.join(",", cells);
  }

  private static String decimal(double value) {
    return Double.isNaN(value) ? "" : String.format(Locale.ROOT, "%.8f", value);
  }

  private static String csv(String value) {
    return '"' + value.replace("\"", "\"\"").replace("\n", " ") + '"';
  }
}
