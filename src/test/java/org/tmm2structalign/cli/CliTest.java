package org.tmm2structalign.cli;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.tmm2structalign.run.BaselineScorer.BaselineScore;
import org.tmm2structalign.testing.SyntheticDataset;

@DisplayName("Unit tests for the command-line interface")
class CliTest {

  @TempDir Path directory;

  private Path dataRoot;
  private Path outputRoot;

  @BeforeEach
  void setUp() throws IOException {
    dataRoot = SyntheticDataset.write(directory.resolve("GPCRdb"));
    outputRoot = directory.resolve("results");
  }

  @Nested
  @DisplayName("When using the run command")
  class Run {

    @Test
    @DisplayName("given a small budget, when running, then results are written and status is 0")
    void givenSmallBudget_whenRunning_thenResultsAreWritten() {
      // Act
      int status = TmStructAlignCli.execute(runArguments("run1", "left,right,true"));

      // Assert
      Path runDirectory = outputRoot.resolve("ejecuciones/synthetic/run1");
      assertEquals(0, status);
      assertTrue(Files.isRegularFile(runDirectory.resolve("FUN.tsv")));
      assertTrue(Files.isRegularFile(runDirectory.resolve("metadata.json")));
    }

    @Test
    @DisplayName("given an existing run directory, when running, then status is 2")
    void givenExistingRunDirectory_whenRunning_thenStatusIsTwo() throws IOException {
      // Arrange
      Files.createDirectories(outputRoot.resolve("ejecuciones/synthetic/run2"));

      // Act & Assert
      assertEquals(2, TmStructAlignCli.execute(runArguments("run2", "left,right")));
    }

    @Test
    @DisplayName("given an unknown seed method, when running, then an error report is written")
    void givenUnknownSeedMethod_whenRunning_thenErrorReportIsWritten() throws IOException {
      // Act
      int status = TmStructAlignCli.execute(runArguments("run3", "left,unknown"));

      // Assert
      Path error = outputRoot.resolve("ejecuciones/synthetic/run3").resolve(RunCommand.ERROR_FILE);
      assertEquals(1, status);
      assertTrue(Files.readString(error).contains("unknown"));
    }

    @Test
    @DisplayName("given missing options, when running, then picocli reports a usage error")
    void givenMissingOptions_whenRunning_thenUsageErrorIsReported() {
      // Act & Assert
      assertEquals(2, TmStructAlignCli.execute("run", "--dataset", "x"));
      assertEquals(0, TmStructAlignCli.execute("--version"));
    }

    private String[] runArguments(String runId, String seeds) {
      return new String[] {
        "run",
        "--data-root", dataRoot.toString(),
        "--dataset", SyntheticDataset.NAME,
        "--seeds", seeds,
        "--output-root", outputRoot.toString(),
        "--run-id", runId,
        "--seed", "7",
        "--max-evaluations", "60",
        "--population-size", "6",
        "--progress-interval", "20"
      };
    }
  }

  @Nested
  @DisplayName("When using the baseline command")
  class Baseline {

    @Test
    @DisplayName("given a dataset, when scoring, then a CSV file is written")
    void givenDataset_whenScoring_thenCsvFileIsWritten() throws IOException {
      // Arrange
      Path csv = directory.resolve("baselines.csv");

      // Act
      int status =
          TmStructAlignCli.execute(
              "baseline",
              "--data-root", dataRoot.toString(),
              "--datasets", SyntheticDataset.NAME,
              "--output", csv.toString());

      // Assert
      List<String> lines = Files.readAllLines(csv);
      assertEquals(0, status);
      assertEquals("dataset,method,SOPwTP,lDDT,file,status,message", lines.get(0));
      assertEquals(4, lines.size());
    }

    @Test
    @DisplayName("given a failed score, when formatting, then status and message are reported")
    void givenFailedScore_whenFormatting_thenStatusAndMessageAreReported() {
      // Arrange
      BaselineScore failed = new BaselineScore("d", "m", "f.fasta", List.of(), "bad \"file\"");

      // Act
      String csv = BaselineCommand.toCsv(List.of(failed));

      // Assert
      assertEquals("\"d\",\"m\",,,\"f.fasta\",failed,\"bad \"\"file\"\"\"", csv.split("\n")[1]);
    }
  }
}
