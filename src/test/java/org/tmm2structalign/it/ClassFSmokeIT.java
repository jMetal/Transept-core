package org.tmm2structalign.it;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assumptions.assumeTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.tmm2structalign.algorithm.AlgorithmConfiguration;
import org.tmm2structalign.objective.LddtAggregation;
import org.tmm2structalign.run.ExperimentRunner;
import org.tmm2structalign.run.ExperimentRunner.RunSummary;
import org.tmm2structalign.run.RunConfiguration;
import org.tmm2structalign.structure.MappingPolicy;

/** Short run on the real classF dataset. Set {@code GPCR_DATASET} to the GPCRdb directory. */
@DisplayName("Integration tests: smoke run on the real classF dataset")
class ClassFSmokeIT {

  @TempDir Path directory;

  @Test
  @DisplayName("given classF and sequence-tool seeds, when running briefly, then it completes")
  void givenClassF_whenRunningBriefly_thenItCompletes() throws IOException {
    // Arrange
    String dataset = System.getenv("GPCR_DATASET");
    assumeTrue(dataset != null && Files.isDirectory(Path.of(dataset)), "GPCR_DATASET not set");
    RunConfiguration configuration =
        new RunConfiguration(
            Path.of(dataset),
            "classF",
            List.of("clustalw", "kalign", "mafft", "muscle5", "tcoffee", "tmaligner"),
            directory.resolve("run"),
            1L,
            AlgorithmConfiguration.withDefaults(200, 10, 2),
            MappingPolicy.VERIFIED,
            LddtAggregation.HARMONIC_OVER_LENGTH,
            15f,
            50);

    // Act
    RunSummary summary = new ExperimentRunner().run(configuration);

    // Assert
    assertEquals(10, summary.frontSize());
    assertTrue(summary.evaluations() >= 200);
    assertTrue(Files.readString(summary.runDirectory().resolve("metadata.json"))
        .contains("\"proteinsWithoutStructure\":[]"));
  }
}
