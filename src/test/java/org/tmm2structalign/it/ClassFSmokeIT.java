package org.tmm2structalign.it;

import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assumptions.assumeTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.io.TempDir;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.EnumSource;
import org.tmm2structalign.algorithm.AlgorithmConfiguration;
import org.tmm2structalign.algorithm.AlgorithmType;
import org.tmm2structalign.objective.LddtAggregation;
import org.tmm2structalign.run.ExperimentRunner;
import org.tmm2structalign.run.ExperimentRunner.RunSummary;
import org.tmm2structalign.run.RunConfiguration;
import org.tmm2structalign.structure.MappingPolicy;

/** Short run on the real classF dataset. Set {@code GPCR_DATASET} to the GPCRdb directory. */
@DisplayName("Integration tests: smoke run on the real classF dataset")
class ClassFSmokeIT {

  @TempDir Path directory;

  @ParameterizedTest
  @EnumSource(AlgorithmType.class)
  @DisplayName("given classF and an algorithm, when running briefly, then it completes")
  void givenClassF_whenRunningBriefly_thenItCompletes(AlgorithmType type) throws IOException {
    // Arrange
    String dataset = System.getenv("GPCR_DATASET");
    assumeTrue(dataset != null && Files.isDirectory(Path.of(dataset)), "GPCR_DATASET not set");
    RunConfiguration configuration =
        new RunConfiguration(
            Path.of(dataset),
            "classF",
            List.of("clustalw", "kalign", "mafft", "muscle5", "tcoffee", "tmaligner"),
            directory.resolve(type.name()),
            1L,
            new AlgorithmConfiguration(
                200, 12, 12, 0.8, type == AlgorithmType.PAES ? 1.0 : 0.2, 2),
            type,
            MappingPolicy.VERIFIED,
            LddtAggregation.HARMONIC_OVER_LENGTH,
            15f,
            50);

    // Act
    RunSummary summary = new ExperimentRunner().run(configuration);

    // Assert
    assertTrue(summary.frontSize() > 0);
    assertTrue(summary.evaluations() >= 200);
    assertTrue(Files.readString(summary.runDirectory().resolve("metadata.json"))
        .contains("\"proteinsWithoutStructure\":[]"));
  }
}
