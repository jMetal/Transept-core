package org.transept.it;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.transept.algorithm.AlgorithmConfiguration;
import org.transept.alignment.GapColumnCompactor;
import org.transept.alignment.Msa;
import org.transept.io.FastaAlignmentReader;
import org.transept.io.ThreeLineTopologyReader;
import org.transept.model.ProteinSet;
import org.transept.run.ExperimentRunner;
import org.transept.run.ExperimentRunner.RunSummary;
import org.transept.run.RunConfiguration;
import org.transept.run.RunDirectoryExistsException;
import org.transept.structure.MappingPolicy;
import org.transept.testing.SyntheticDataset;

@DisplayName("Integration tests: full runs on a synthetic dataset")
class SyntheticEndToEndIT {

  private static final int POPULATION = 10;

  @TempDir Path directory;

  private Path dataRoot;
  private ExperimentRunner runner;

  @BeforeEach
  void setUp() throws IOException {
    dataRoot = SyntheticDataset.write(directory.resolve("GPCRdb"));
    runner = new ExperimentRunner();
  }

  @Test
  @DisplayName("given a synthetic dataset, when running, then every output is consistent")
  void givenSyntheticDataset_whenRunning_thenEveryOutputIsConsistent() throws IOException {
    // Arrange
    RunConfiguration configuration = configuration("run", 500);

    // Act
    RunSummary summary = runner.run(configuration);

    // Assert
    Path run = summary.runDirectory();
    List<String> fun = Files.readAllLines(run.resolve("FUN.csv"));
    assertEquals(POPULATION, summary.frontSize());
    assertEquals(POPULATION, fun.size());
    assertTrue(summary.evaluations() >= 500);
    ProteinSet proteins = new ThreeLineTopologyReader().read(configuration.topologyFile());
    for (int index = 0; index < fun.size(); index++) {
      assertEquals(2, fun.get(index).split(",").length);
      Msa msa = new FastaAlignmentReader().read(run.resolve("MSASol" + index + ".fasta"), proteins);
      assertEquals(0, GapColumnCompactor.countGapOnlyColumns(msa.rows()), "MSASol" + index);
    }
    for (String file : List.of("config.json", "metadata.json", "progress.jsonl", "seed.txt",
        "runtime.csv")) {
      assertTrue(Files.isRegularFile(run.resolve(file)), file);
    }
    List<String> runtime = Files.readAllLines(run.resolve("runtime.csv"));
    assertEquals("runtime_ms,runtime_seconds,evaluations", runtime.get(0));
    assertEquals(String.valueOf(summary.evaluations()), runtime.get(1).split(",")[2]);
    String metadata = Files.readString(run.resolve("metadata.json"));
    assertTrue(metadata.contains("\"proteinsWithoutStructure\":[]"));
    assertTrue(metadata.contains("\"status\":\"NUMBERED\""));
    assertTrue(Files.readAllLines(run.resolve("progress.jsonl")).size() >= 4);
  }

  @Test
  @DisplayName("given an existing run, when running again, then nothing is overwritten")
  void givenExistingRun_whenRunningAgain_thenNothingIsOverwritten() throws IOException {
    // Arrange
    RunConfiguration configuration = configuration("again", 40);
    runner.run(configuration);
    String fun = Files.readString(configuration.runDirectory().resolve("FUN.csv"));

    // Act & Assert
    assertThrows(RunDirectoryExistsException.class, () -> runner.run(configuration));
    assertEquals(fun, Files.readString(configuration.runDirectory().resolve("FUN.csv")));
  }

  private RunConfiguration configuration(String runId, int evaluations) {
    return new RunConfiguration(
        dataRoot,
        SyntheticDataset.NAME,
        SyntheticDataset.SEED_METHODS,
        directory.resolve("out").resolve(runId),
        3L,
        AlgorithmConfiguration.withDefaults(evaluations, POPULATION, 1),
        MappingPolicy.VERIFIED,
        15f,
        100);
  }
}
