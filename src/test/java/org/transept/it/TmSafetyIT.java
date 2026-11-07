package org.transept.it;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.io.TempDir;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.EnumSource;
import org.transept.algorithm.AlgorithmConfiguration;
import org.transept.algorithm.AlgorithmType;
import org.transept.algorithm.MutationType;
import org.transept.alignment.Msa;
import org.transept.io.FastaAlignmentReader;
import org.transept.io.ThreeLineTopologyReader;
import org.transept.model.ProteinSet;
import org.transept.objective.ObjectiveSet;
import org.transept.run.ExperimentRunner;
import org.transept.run.RunConfiguration;
import org.transept.structure.MappingPolicy;
import org.transept.testing.Alignments;
import org.transept.testing.SyntheticDataset;

/**
 * With seeds that have no gaps inside helices, the TM-safe crossover and the TM-safe mutation
 * must never produce an alignment with a gap run between two transmembrane residues.
 */
@DisplayName("Integration tests: transmembrane helices are never broken with TM-safe operators")
class TmSafetyIT {

  @TempDir Path directory;

  private Path dataRoot;

  @BeforeEach
  void setUp() throws IOException {
    dataRoot = SyntheticDataset.write(directory.resolve("GPCRdb"));
  }

  @ParameterizedTest
  @EnumSource(
      value = AlgorithmType.class,
      names = {"RANDOM_SEARCH"},
      mode = EnumSource.Mode.EXCLUDE)
  @DisplayName("given intact seeds, when running, then no final alignment breaks a helix")
  void givenIntactSeeds_whenRunning_thenNoFinalAlignmentBreaksHelix(AlgorithmType type)
      throws IOException {
    // Arrange: "left" and "right" only have gaps at the ends of the sequences
    Path run = directory.resolve("out").resolve(type.name());
    RunConfiguration configuration =
        new RunConfiguration(
            dataRoot,
            SyntheticDataset.NAME,
            List.of("left", "right"),
            run,
            7L,
            new AlgorithmConfiguration(
                400, 8, 8, 0.9, 1.0, 1, MutationType.TM_SAFE_GAP_SHIFT),
            type,
            MappingPolicy.VERIFIED,
            15f,
            100,
            type.isSingleObjective()
                ? ObjectiveSet.LDDT
                : ObjectiveSet.SOPWTP_LDDT);

    // Act
    int solutions = new ExperimentRunner().run(configuration).frontSize();

    // Assert
    assertEquals(0, breakingRuns(configuration, run, solutions), type + " broke a helix");
    assertEquals(
        true, Files.readString(run.resolve("config.json")).contains("TM_SAFE_GAP_SHIFT"));
  }

  @ParameterizedTest
  @EnumSource(MutationType.class)
  @DisplayName("given a mutation, when running, then no helix breaks")
  void givenMutation_whenRunning_thenNoHelixBreaks(MutationType mutation) throws IOException {
    // Arrange
    Path run = directory.resolve("operators").resolve(mutation.toString());
    RunConfiguration configuration =
        new RunConfiguration(
            dataRoot,
            SyntheticDataset.NAME,
            List.of("left", "right"),
            run,
            11L,
            new AlgorithmConfiguration(
                400, 8, 8, 0.9, 1.0, 1, mutation),
            AlgorithmType.NSGA_II,
            MappingPolicy.VERIFIED,
            15f,
            100);

    // Act
    int solutions = new ExperimentRunner().run(configuration).frontSize();

    // Assert
    assertEquals(
        0, breakingRuns(configuration, run, solutions), mutation.toString());
  }

  @ParameterizedTest
  @EnumSource(
      value = AlgorithmType.class,
      names = {"NSGA_II", "SPEA2", "GENETIC_ALGORITHM", "SMS_EMOA", "MOEA_D", "PAES"})
  @DisplayName("given a limited shift and several cuts, when running, then no helix breaks")
  void givenLimitedShiftAndSeveralCuts_whenRunning_thenNoHelixBreaks(AlgorithmType type)
      throws IOException {
    // Arrange
    Path run = directory.resolve("settings").resolve(type.name());
    RunConfiguration configuration =
        new RunConfiguration(
            dataRoot,
            SyntheticDataset.NAME,
            List.of("left", "right"),
            run,
            13L,
            new AlgorithmConfiguration(
                400, 8, 8, 0.9, 1.0, 1, MutationType.TM_SAFE_GAP_SHIFT, null, "", 3, 4),
            type,
            MappingPolicy.VERIFIED,
            15f,
            100,
            type.isSingleObjective() ? ObjectiveSet.LDDT : ObjectiveSet.SOPWTP_LDDT);

    // Act
    int solutions = new ExperimentRunner().run(configuration).frontSize();

    // Assert
    assertEquals(0, breakingRuns(configuration, run, solutions), type + " broke a helix");
  }

  private static int breakingRuns(RunConfiguration configuration, Path run, int solutions)
      throws IOException {
    ProteinSet proteins = new ThreeLineTopologyReader().read(configuration.topologyFile());
    int breaking = 0;
    for (int index = 0; index < solutions; index++) {
      Msa msa =
          new FastaAlignmentReader().read(run.resolve("MSASol" + index + ".fasta"), proteins);
      for (int sequence = 0; sequence < msa.size(); sequence++) {
        breaking += Alignments.helixBreakingRuns(proteins.get(sequence), msa.rows().get(sequence));
      }
    }
    return breaking;
  }
}
