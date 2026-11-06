package org.transept.evaluation;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class ScoreReportTest {
  @TempDir
  Path tempDirectory;

  @Test
  void recursivelyDiscoversSeedAlignmentFilesInStableOrder() throws Exception {
    Path seed2 = Files.createDirectories(tempDirectory.resolve("seed_002"));
    Path seed1 = Files.createDirectories(tempDirectory.resolve("seed_001"));
    Files.writeString(seed2.resolve("MSA2.fasta"), ">a\nA\n");
    Files.writeString(seed1.resolve("MSA1.fasta"), ">a\nA\n");
    Files.writeString(seed1.resolve("notes.txt"), "ignored");

    List<Path> result = ScoreReport.discoverAlignments(tempDirectory, "MSA*.fasta");

    assertEquals(List.of(
        seed1.resolve("MSA1.fasta"),
        seed2.resolve("MSA2.fasta")), result);
  }

  @Test
  void reportsAnAlignmentThatCannotBeComparedAndGoesOn() throws Exception {
    Path reference = tempDirectory.resolve("reference.fasta");
    Path topology = tempDirectory.resolve("topology.3line");
    Path root = Files.createDirectories(tempDirectory.resolve("alignments"));
    Path output = tempDirectory.resolve("report.csv");
    Files.writeString(reference, ">a\nACD\n>b\nACD\n");
    Files.writeString(topology, ">a\nACD\nMMO\n>b\nACD\nMMO\n");
    Files.writeString(root.resolve("MSA0.fasta"), ">a\nACD\n>b\nACD\n");
    Files.writeString(root.resolve("MSA1.fasta"), ">a\nACD\n>other\nACD\n");

    ScoreReport.Summary summary =
        ScoreReport.write(reference, topology, root, "MSA*.fasta", output);

    assertEquals(new ScoreReport.Summary(1, 1), summary);
    List<String> lines = Files.readAllLines(output);
    assertEquals(3, lines.size());
    int columns = lines.get(0).split(",").length;
    assertEquals(columns, lines.get(1).split(",", -1).length);
    assertTrue(lines.get(1).contains(",ok,"));
    assertTrue(lines.get(2).contains(",failed,"));
    assertEquals(columns, csvColumns(lines.get(2)));
  }

  @Test
  void refusesToOverwriteAnInput() throws Exception {
    Path reference = tempDirectory.resolve("reference.fasta");
    Path topology = tempDirectory.resolve("topology.3line");
    Files.writeString(reference, ">a\nA\n");
    Files.writeString(topology, ">a\nA\nM\n");
    Path alignment = tempDirectory.resolve("MSA0.fasta");
    Files.writeString(alignment, ">a\nA\n");

    assertThrows(
        IOException.class,
        () -> ScoreReport.write(reference, topology, tempDirectory, "MSA*.fasta", alignment));
  }

  @Test
  void failsWhenNoAlignmentMatches() throws Exception {
    Path file = tempDirectory.resolve("reference.fasta");
    Files.writeString(file, ">a\nA\n");

    assertThrows(
        IOException.class,
        () ->
            ScoreReport.write(
                file, file, tempDirectory, "MSA*.fasta", tempDirectory.resolve("o")));
  }

  private static int csvColumns(String line) {
    int columns = 1;
    boolean quoted = false;
    for (char character : line.toCharArray()) {
      if (character == '"') {
        quoted = !quoted;
      } else if (character == ',' && !quoted) {
        columns++;
      }
    }
    return columns;
  }
}
