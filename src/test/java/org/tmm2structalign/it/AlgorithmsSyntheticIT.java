package org.tmm2structalign.it;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.io.TempDir;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.EnumSource;
import org.tmm2structalign.algorithm.AlgorithmConfiguration;
import org.tmm2structalign.algorithm.AlgorithmType;
import org.tmm2structalign.alignment.GapColumnCompactor;
import org.tmm2structalign.alignment.Msa;
import org.tmm2structalign.io.FastaAlignmentReader;
import org.tmm2structalign.io.ThreeLineTopologyReader;
import org.tmm2structalign.model.ProteinSet;
import org.tmm2structalign.objective.LddtAggregation;
import org.tmm2structalign.run.ExperimentRunner;
import org.tmm2structalign.run.ExperimentRunner.RunSummary;
import org.tmm2structalign.run.RunConfiguration;
import org.tmm2structalign.structure.MappingPolicy;
import org.tmm2structalign.testing.SyntheticDataset;

@DisplayName("Integration tests: every algorithm on the synthetic dataset")
class AlgorithmsSyntheticIT {

  private static final int POPULATION = 8;
  private static final int EVALUATIONS = 120;

  @TempDir Path directory;

  private Path dataRoot;

  @BeforeEach
  void setUp() throws IOException {
    dataRoot = SyntheticDataset.write(directory.resolve("GPCRdb"));
  }

  @ParameterizedTest
  @EnumSource(AlgorithmType.class)
  @DisplayName("given an algorithm, when running, then the outputs are complete and valid")
  void givenAlgorithm_whenRunning_thenOutputsAreCompleteAndValid(AlgorithmType type)
      throws IOException {
    // Act
    RunSummary summary = run(type, "run", 4L);

    // Assert
    Path run = summary.runDirectory();
    List<String> fun = Files.readAllLines(run.resolve("FUN.csv"));
    assertFalse(fun.isEmpty());
    assertEquals(summary.frontSize(), fun.size());
    assertTrue(summary.evaluations() >= EVALUATIONS);
    boolean unboundedFront =
        type == AlgorithmType.RANDOM_SEARCH || type == AlgorithmType.RANDOM_SEEDED;
    assertTrue(unboundedFront || fun.size() <= POPULATION, type + ": " + fun.size());
    ProteinSet proteins =
        new ThreeLineTopologyReader()
            .read(dataRoot.resolve("sequences/tmregions/synthetic_predicted_topologies.3line"));
    for (int index = 0; index < fun.size(); index++) {
      Msa msa =
          new FastaAlignmentReader().read(run.resolve("MSASol" + index + ".fasta"), proteins);
      assertEquals(0, GapColumnCompactor.countGapOnlyColumns(msa.rows()), type + " " + index);
    }
    assertTrue(Files.readString(run.resolve("metadata.json"))
        .contains("\"algorithmType\":\"" + type + "\""));
  }

  @ParameterizedTest
  @EnumSource(AlgorithmType.class)
  @DisplayName("given an algorithm and a seed, when running twice, then the results are identical")
  void givenAlgorithmAndSeed_whenRunningTwice_thenResultsAreIdentical(AlgorithmType type)
      throws IOException {
    // Act
    Path first = run(type, "first", 9L).runDirectory();
    Path second = run(type, "second", 9L).runDirectory();

    // Assert
    assertEquals(
        Files.readString(first.resolve("FUN.csv")), Files.readString(second.resolve("FUN.csv")));
    assertEquals(
        Files.readString(first.resolve("MSASol0.fasta")),
        Files.readString(second.resolve("MSASol0.fasta")));
  }

  private RunSummary run(AlgorithmType type, String runId, long seed) throws IOException {
    return new ExperimentRunner()
        .run(
            new RunConfiguration(
                dataRoot,
                SyntheticDataset.NAME,
                List.of("left", "right"),
                directory.resolve("out").resolve(type.name()).resolve(runId),
                seed,
                new AlgorithmConfiguration(
                    EVALUATIONS,
                    POPULATION,
                    POPULATION,
                    0.8,
                    type == AlgorithmType.PAES ? 1.0 : 0.2,
                    1),
                type,
                MappingPolicy.VERIFIED,
                LddtAggregation.HARMONIC_OVER_LENGTH,
                15f,
                40));
  }
}
