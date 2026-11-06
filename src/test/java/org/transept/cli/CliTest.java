package org.transept.cli;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.transept.algorithm.AlgorithmType;
import org.transept.algorithm.MutationType;
import org.transept.run.BaselineScorer.BaselineScore;
import org.transept.testing.SyntheticDataset;

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
  @DisplayName("When using the score command")
  class Score {

    private Path write(String name, String content) throws IOException {
      Path path = directory.resolve(name);
      Files.createDirectories(path.getParent());
      Files.writeString(path, content);
      return path;
    }

    @Test
    @DisplayName("given alignments and a reference, when scoring, then a report is written")
    void givenAlignments_whenScoring_thenReportIsWritten() throws IOException {
      // Arrange
      Path reference = write("reference.fasta", ">a\nACD\n>b\nACD\n");
      Path topology = write("topology.3line", ">a\nACD\nMMO\n>b\nACD\nMMO\n");
      write("alignments/mafft_x.fasta", ">a\nACD\n>b\nACD\n");
      Path report = directory.resolve("report.csv");

      // Act
      int status =
          TranseptCli.execute(
              "score",
              "--reference", reference.toString(),
              "--topology", topology.toString(),
              "--msa-root", directory.resolve("alignments").toString(),
              "--glob", "*.fasta",
              "--output", report.toString());

      // Assert
      assertEquals(0, status);
      assertTrue(Files.readString(report).contains("mafft_x.fasta"));
    }

    @Test
    @DisplayName("given no alignment can be compared, when scoring, then the status is one")
    void givenOnlyWrongAlignments_whenScoring_thenStatusIsOne() throws IOException {
      // Arrange
      Path reference = write("reference.fasta", ">a\nACD\n>b\nACD\n");
      Path topology = write("topology.3line", ">a\nACD\nMMO\n>b\nACD\nMMO\n");
      write("alignments/MSA0.fasta", ">a\nACD\n>other\nACD\n");

      // Act
      int status =
          TranseptCli.execute(
              "score",
              "--reference", reference.toString(),
              "--topology", topology.toString(),
              "--msa-root", directory.resolve("alignments").toString(),
              "--output", directory.resolve("report.csv").toString());

      // Assert
      assertEquals(1, status);
    }

    @Test
    @DisplayName("given a missing reference, when scoring, then a usage error is reported")
    void givenMissingReference_whenScoring_thenUsageErrorIsReported() throws IOException {
      // Arrange
      Path topology = write("topology.3line", ">a\nACD\nMMO\n");

      // Act
      int status =
          TranseptCli.execute(
              "score",
              "--reference", directory.resolve("gone.fasta").toString(),
              "--topology", topology.toString(),
              "--msa-root", directory.toString());

      // Assert
      assertEquals(2, status);
    }
  }

  @Nested
  @DisplayName("When using the run command")
  class Run {

    @Test
    @DisplayName("given a small budget, when running, then results are written and status is 0")
    void givenSmallBudget_whenRunning_thenResultsAreWritten() {
      // Act
      int status = TranseptCli.execute(runArguments("run1", "left,right,true"));

      // Assert
      Path runDirectory = outputRoot.resolve("ejecuciones/synthetic/run1");
      assertEquals(0, status);
      assertTrue(Files.isRegularFile(runDirectory.resolve("FUN.csv")));
      assertTrue(Files.isRegularFile(runDirectory.resolve("metadata.json")));
    }

    @Test
    @DisplayName("given an existing run directory, when running, then status is 2")
    void givenExistingRunDirectory_whenRunning_thenStatusIsTwo() throws IOException {
      // Arrange
      Files.createDirectories(outputRoot.resolve("ejecuciones/synthetic/run2"));

      // Act & Assert
      assertEquals(2, TranseptCli.execute(runArguments("run2", "left,right")));
    }

    @Test
    @DisplayName("given an unknown seed method, when running, then an error report is written")
    void givenUnknownSeedMethod_whenRunning_thenErrorReportIsWritten() throws IOException {
      // Act
      int status = TranseptCli.execute(runArguments("run3", "left,unknown"));

      // Assert
      Path error = outputRoot.resolve("ejecuciones/synthetic/run3").resolve(RunCommand.ERROR_FILE);
      assertEquals(1, status);
      assertTrue(Files.readString(error).contains("unknown"));
    }

    @Test
    @DisplayName("given missing options, when running, then picocli reports a usage error")
    void givenMissingOptions_whenRunning_thenUsageErrorIsReported() {
      // Act & Assert
      assertEquals(2, TranseptCli.execute("run", "--dataset", "x"));
      assertEquals(0, TranseptCli.execute("--version"));
    }

    @Test
    @DisplayName("given internal and repaired seeds, when running, then the seeds are recorded")
    void givenInternalAndRepairedSeeds_whenRunning_thenTheSeedsAreRecorded() throws IOException {
      // Act: the synthetic seeds left and right are repaired too
      int status =
          TranseptCli.execute(
              "run",
              "--data-root", dataRoot.toString(),
              "--dataset", SyntheticDataset.NAME,
              "--seeds", "internal,left",
              "--repair-seeds",
              "--output-root", outputRoot.toString(),
              "--run-id", "internal",
              "--seed", "3",
              "--max-evaluations", "50",
              "--population-size", "8");

      // Assert
      Path run = outputRoot.resolve("ejecuciones/synthetic/internal");
      String metadata = Files.readString(run.resolve("metadata.json"));
      assertEquals(0, status);
      assertTrue(metadata.contains("\"method\":\"internal:biojava-g20x2\""), metadata);
      assertTrue(metadata.contains("\"repairedGapRuns\":"), metadata);
      assertTrue(Files.isRegularFile(run.resolve("seeds/biojava-g20x2.fasta")));
      assertTrue(Files.readString(run.resolve("config.json")).contains("\"repairSeeds\":true"));
    }

    @Test
    @DisplayName("given seeds computed in advance, when running, then the run reads them")
    void givenSeedsComputedInAdvance_whenRunning_thenTheRunReadsThem() throws IOException {
      // Arrange
      int prepared =
          TranseptCli.execute(
              "seeds",
              "--data-root", dataRoot.toString(),
              "--datasets", SyntheticDataset.NAME,
              "--output-root", outputRoot.toString(),
              "--count", "3");

      // Act: the synthetic proteins are so short that every setting gives the same alignment, so
      // a seed of the dataset is added to fill the population
      int status =
          TranseptCli.execute(
              "run",
              "--data-root", dataRoot.toString(),
              "--dataset", SyntheticDataset.NAME,
              "--output-root", outputRoot.toString(),
              "--seeds", "internal,left",
              "--run-id", "precomputed",
              "--seed", "3",
              "--max-evaluations", "50",
              "--population-size", "8");

      // Assert
      assertEquals(0, prepared);
      assertEquals(0, status);
      assertTrue(
          Files.isRegularFile(outputRoot.resolve("internal_seeds/synthetic/manifest.tsv")));
      Path run = outputRoot.resolve("ejecuciones/synthetic/precomputed");
      String metadata = Files.readString(run.resolve("metadata.json"));
      assertTrue(metadata.contains("\"origin\":\"precomputed\""), metadata);
      assertTrue(!metadata.contains("\"origin\":\"computed\""), metadata);
      assertEquals(
          1, TranseptCli.execute(
              "seeds", "--data-root", dataRoot.toString(), "--datasets", "missing",
              "--output-root", outputRoot.toString()));
    }

    @Test
    @DisplayName("given a time limit, when running, then it stops by time and records it")
    void givenTimeLimit_whenRunning_thenItStopsByTimeAndRecordsIt() throws IOException {
      // Act
      int status =
          TranseptCli.execute(
              "run",
              "--data-root", dataRoot.toString(),
              "--dataset", SyntheticDataset.NAME,
              "--seeds", "left,right",
              "--output-root", outputRoot.toString(),
              "--run-id", "timed",
              "--seed", "3",
              "--max-time", "0.5",
              "--population-size", "6");

      // Assert
      Path run = outputRoot.resolve("ejecuciones/synthetic/timed");
      String config = Files.readString(run.resolve("config.json"));
      assertEquals(0, status);
      assertTrue(config.contains("\"maxTimeSeconds\":0.5"), config);
      assertTrue(config.contains("\"maxEvaluations\":null"), config);
    }

    @Test
    @DisplayName("given both stopping conditions, when running, then it fails with status 1")
    void givenBothStoppingConditions_whenRunning_thenStatusIsOne() {
      // Arrange
      String[] arguments = runArguments("both", "left,right");
      String[] withTime = Arrays.copyOf(arguments, arguments.length + 2);
      withTime[arguments.length] = "--max-time";
      withTime[arguments.length + 1] = "5";

      // Act & Assert
      assertEquals(1, TranseptCli.execute(withTime));
    }

    @Test
    @DisplayName("given a non-positive time, when resolving the limit, then it fails")
    void givenNonPositiveTime_whenResolvingLimit_thenIllegalArgumentExceptionIsThrown() {
      // Arrange
      RunCommand command = new RunCommand();
      command.maxTime = 0.0;

      // Act & Assert
      assertThrows(IllegalArgumentException.class, command::maxRunTime);
      command.maxTime = null;
      assertNull(command.maxRunTime());
    }

    @Test
    @DisplayName("given a configuration, when running NSGA-II, then metadata records it")
    void givenConfiguration_whenRunningNsgaTwo_thenMetadataRecordsIt() throws IOException {
      // Arrange
      Path file = directory.resolve("configuration.txt");
      Files.writeString(file, "\n--selection random --mutationProbability 0.5\n");
      String[] arguments = runArguments("configured", "left,right");
      String[] withFile = Arrays.copyOf(arguments, arguments.length + 2);
      withFile[arguments.length] = "--configuration-file";
      withFile[arguments.length + 1] = file.toString();

      // Act
      int status = TranseptCli.execute(withFile);

      // Assert
      String metadata =
          Files.readString(outputRoot.resolve("ejecuciones/synthetic/configured/metadata.json"));
      assertEquals(0, status);
      assertTrue(metadata.contains("--selection random"), metadata);
      assertTrue(metadata.contains("--mutationProbability 0.5"), metadata);
    }

    @Test
    @DisplayName("given an option and a configuration setting the same parameter, then it fails")
    void givenOptionAndConfigurationSettingSameParameter_whenResolving_thenItFails() {
      // Arrange
      RunCommand command = new RunCommand();
      command.mutation = MutationType.MERGE_GAP_RUNS;
      command.configuration = "--mutation randomChoice";

      // Act & Assert
      assertThrows(IllegalArgumentException.class, command::configurationText);
      command.mutation = null;
      command.configurationFile = directory.resolve("x.txt");
      assertThrows(IllegalArgumentException.class, command::configurationText);
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
  @DisplayName("When using the algorithms command")
  class Algorithms {

    @Test
    @DisplayName("given the command, when run, then it lists every algorithm with its objectives")
    void givenCommand_whenRun_thenItListsEveryAlgorithmWithItsObjectives() {
      // Arrange
      PrintStream original = System.out;
      ByteArrayOutputStream output = new ByteArrayOutputStream();
      System.setOut(new PrintStream(output, true, StandardCharsets.UTF_8));

      // Act
      int status;
      try {
        status = TranseptCli.execute("algorithms");
      } finally {
        System.setOut(original);
      }

      // Assert
      String json = output.toString(StandardCharsets.UTF_8);
      assertEquals(0, status);
      for (AlgorithmType type : AlgorithmType.values()) {
        assertTrue(json.contains("\"name\":\"" + type + "\""), type.toString());
      }
      assertTrue(
          json.contains(
              "\"name\":\"GENETIC_ALGORITHM\",\"minObjectives\":1,\"maxObjectives\":1"),
          json);
      assertTrue(json.contains("\"parameterSpace\":\"MSASPEA2.yaml\""), json);
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
          TranseptCli.execute(
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

  @Nested
  @DisplayName("When choosing the algorithm")
  class Algorithm {

    @Test
    @DisplayName("given PAES, when running, then the metadata records it and pm defaults to 1")
    void givenPaes_whenRunning_thenMetadataRecordsItAndMutationIsOne() throws IOException {
      // Act
      int status =
          TranseptCli.execute(
              "run",
              "--data-root", dataRoot.toString(),
              "--dataset", SyntheticDataset.NAME,
              "--seeds", "left,right",
              "--output-root", outputRoot.toString(),
              "--run-id", "paes",
              "--seed", "3",
              "--max-evaluations", "50",
              "--population-size", "8",
              "--algorithm", "PAES");

      // Assert
      Path runDirectory = outputRoot.resolve("ejecuciones/synthetic/paes");
      assertEquals(0, status);
      assertTrue(Files.readString(runDirectory.resolve("metadata.json"))
          .contains("\"algorithmType\":\"PAES\""));
      assertTrue(Files.readString(runDirectory.resolve("config.json"))
          .contains("\"mutationProbability\":1.0"));
    }

    @Test
    @DisplayName("given new operators, when running, then the configuration records them")
    void givenNewOperators_whenRunning_thenConfigurationRecordsThem() throws IOException {
      // Act
      int status =
          TranseptCli.execute(
              "run",
              "--data-root", dataRoot.toString(),
              "--dataset", SyntheticDataset.NAME,
              "--seeds", "left,right",
              "--output-root", outputRoot.toString(),
              "--run-id", "operators",
              "--seed", "3",
              "--max-evaluations", "50",
              "--population-size", "8",
              "--mutation", "RANDOM_CHOICE");

      // Assert
      String config =
          Files.readString(outputRoot.resolve("ejecuciones/synthetic/operators/config.json"));
      assertEquals(0, status);
      assertTrue(config.contains("\"mutationType\":\"RANDOM_CHOICE\""), config);
    }

    @ParameterizedTest
    @CsvSource({"PAES, , 1.0", "NSGA_II, , 0.2", "SMS_EMOA, , 0.2", "PAES, 0.5, 0.5"})
    @DisplayName("given an algorithm, when resolving pm, then the explicit or default is used")
    void givenAlgorithm_whenResolvingMutationProbability_thenExplicitOrDefaultValueIsUsed(
        AlgorithmType algorithm, Double explicit, double expected) {
      // Arrange
      RunCommand command = new RunCommand();
      command.algorithm = algorithm;
      command.mutationProbability = explicit;

      // Act & Assert
      assertEquals(expected, command.effectiveMutationProbability());
    }

    @Test
    @DisplayName("given an unknown algorithm, when running, then picocli reports a usage error")
    void givenUnknownAlgorithm_whenRunning_thenUsageErrorIsReported() {
      // Act & Assert
      assertEquals(2, TranseptCli.execute("run", "--algorithm", "SPEA2"));
    }
  }
}
