package org.tmm2structalign.it;

import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.tmm2structalign.algorithm.AlgorithmConfiguration;
import org.tmm2structalign.algorithm.AlgorithmType;
import org.tmm2structalign.run.ExperimentRunner;
import org.tmm2structalign.run.RunConfiguration;
import org.tmm2structalign.structure.MappingPolicy;
import org.tmm2structalign.testing.SyntheticDataset;

/**
 * Sanity check: with the same budget, NSGA-II must find better alignments than a random search
 * that places gaps at random.
 */
@DisplayName("Integration tests: NSGA-II against random search")
class SanityCheckIT {

  @TempDir Path directory;

  private Path dataRoot;

  @BeforeEach
  void setUp() throws IOException {
    dataRoot = SyntheticDataset.write(directory.resolve("GPCRdb"));
  }

  @Test
  @DisplayName("given the same budget, when comparing, then NSGA-II beats random search on SOPwTP")
  void givenSameBudget_whenComparing_thenNsgaTwoBeatsRandomSearch() throws IOException {
    // Act
    double nsgaTwo = bestSumOfPairs(AlgorithmType.NSGA_II);
    double randomSearch = bestSumOfPairs(AlgorithmType.RANDOM_SEARCH);

    // Assert
    assertTrue(nsgaTwo > randomSearch, nsgaTwo + " <= " + randomSearch);
  }

  private double bestSumOfPairs(AlgorithmType type) throws IOException {
    Path runDirectory = directory.resolve("out").resolve(type.name());
    new ExperimentRunner()
        .run(
            new RunConfiguration(
                dataRoot,
                SyntheticDataset.NAME,
                List.of("left", "right"),
                runDirectory,
                21L,
                AlgorithmConfiguration.withDefaults(300, 10, 1),
                type,
                MappingPolicy.VERIFIED,
                15f,
                100));
    return Files.readAllLines(runDirectory.resolve("FUN.csv")).stream()
        .mapToDouble(line -> Double.parseDouble(line.split(",")[0]))
        .max()
        .orElseThrow();
  }
}
