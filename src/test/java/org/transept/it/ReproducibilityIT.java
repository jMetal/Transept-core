package org.transept.it;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.transept.algorithm.AlgorithmConfiguration;
import org.transept.run.ExperimentRunner;
import org.transept.run.RunConfiguration;
import org.transept.structure.MappingPolicy;
import org.transept.testing.SyntheticDataset;

/**
 * Seeds exclude the "true" alignment, which dominates every other solution of the synthetic
 * dataset and would make all runs converge to the same population.
 */
@DisplayName("Integration tests: reproducibility from the random seed")
class ReproducibilityIT {

  @TempDir Path directory;

  private Path dataRoot;

  @BeforeEach
  void setUp() throws IOException {
    dataRoot = SyntheticDataset.write(directory.resolve("GPCRdb"));
  }

  @Test
  @DisplayName("given the same seed, when running twice, then results are identical")
  void givenSameSeed_whenRunningTwice_thenResultsAreIdentical() throws IOException {
    // Act
    String first = run("a", 11L, 1);
    String second = run("b", 11L, 1);

    // Assert
    assertEquals(first, second);
  }

  @Test
  @DisplayName("given different seeds, when running, then results differ")
  void givenDifferentSeeds_whenRunning_thenResultsDiffer() throws IOException {
    // Act
    String first = run("c", 11L, 1);
    String second = run("d", 12L, 1);

    // Assert
    assertNotEquals(first, second);
  }

  @Test
  @DisplayName("given one or two threads, when running with the same seed, then results match")
  void givenOneOrTwoThreads_whenRunningWithSameSeed_thenResultsMatch() throws IOException {
    // Act
    String sequential = run("e", 13L, 1);
    String parallel = run("f", 13L, 2);

    // Assert
    assertEquals(sequential, parallel);
  }

  /** Runs and returns FUN.csv followed by every MSASol file, as one string. */
  private String run(String runId, long seed, int threads) throws IOException {
    Path runDirectory = directory.resolve("out").resolve(runId);
    new ExperimentRunner()
        .run(
            new RunConfiguration(
                dataRoot,
                SyntheticDataset.NAME,
                List.of("left", "right"),
                runDirectory,
                seed,
                AlgorithmConfiguration.withDefaults(300, 10, threads),
                MappingPolicy.VERIFIED,
                15f,
                100));
    StringBuilder results = new StringBuilder(Files.readString(runDirectory.resolve("FUN.csv")));
    for (int index = 0; index < 10; index++) {
      results.append(Files.readString(runDirectory.resolve("MSASol" + index + ".fasta")));
    }
    return results.toString();
  }
}
