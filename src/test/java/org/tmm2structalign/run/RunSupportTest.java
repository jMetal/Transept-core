package org.tmm2structalign.run;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.tmm2structalign.algorithm.AlgorithmConfiguration;
import org.tmm2structalign.alignment.MsaSolution;
import org.tmm2structalign.model.ProteinSet;
import org.tmm2structalign.objective.TopologyAwareSumOfPairs;
import org.tmm2structalign.problem.StructuralMsaProblem;
import org.tmm2structalign.structure.MappingPolicy;
import org.tmm2structalign.testing.Alignments;
import org.tmm2structalign.testing.SyntheticDataset;

@DisplayName("Unit tests for Json, RunConfiguration, PrecomputedAlignments, ProgressObserver")
class RunSupportTest {

  @TempDir Path directory;

  @Nested
  @DisplayName("When writing JSON")
  class JsonOutput {

    @Test
    @DisplayName("given nested values, when writing, then valid compact JSON is produced")
    void givenNestedValues_whenWriting_thenCompactJsonIsProduced() {
      // Arrange
      Map<String, Object> value = new LinkedHashMap<>();
      value.put("text", "a\"b\\c\nd\te\u0001");
      value.put("integer", 3);
      value.put("long", 4L);
      value.put("double", 0.5);
      value.put("whole", 2.0);
      value.put("flag", true);
      value.put("none", null);
      value.put("policy", MappingPolicy.VERIFIED);
      value.put("list", List.of(1, "x"));

      // Act
      String json = Json.write(value);

      // Assert
      assertEquals(
          "{\"text\":\"a\\\"b\\\\c\\nd\\te\\u0001\",\"integer\":3,\"long\":4,\"double\":0.5,"
              + "\"whole\":2.0,\"flag\":true,\"none\":null,\"policy\":\"VERIFIED\","
              + "\"list\":[1,\"x\"]}",
          json);
    }

    @Test
    @DisplayName("given unsupported values, when writing, then an exception is thrown")
    void givenUnsupportedValues_whenWriting_thenIllegalArgumentExceptionIsThrown() {
      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> Json.write(new Object()));
      assertThrows(IllegalArgumentException.class, () -> Json.write(Double.NaN));
    }
  }

  @Nested
  @DisplayName("When configuring a run")
  class Configuration {

    @Test
    @DisplayName("given a valid configuration, when deriving paths, then dataset layout is used")
    void givenValidConfiguration_whenDerivingPaths_thenDatasetLayoutIsUsed() {
      // Act
      RunConfiguration configuration = configuration(List.of("mafft"), 100);

      // Assert
      assertEquals(
          directory.resolve("sequences/tmregions/classF_predicted_topologies.3line"),
          configuration.topologyFile());
      assertEquals(directory.resolve("distances/classF"), configuration.distanceDirectory());
      assertEquals(directory.resolve("alphafold_pdb/classF"), configuration.pdbDirectory());
      assertEquals(directory.resolve("precomputed/classF"), configuration.precomputedDirectory());
      assertEquals(List.of("mafft"), configuration.toMap().get("seedMethods"));
      assertEquals("\"VERIFIED\"", Json.write(configuration.toMap().get("mappingPolicy")));
    }

    @Test
    @DisplayName("given invalid values, when creating, then an exception is thrown")
    void givenInvalidValues_whenCreating_thenIllegalArgumentExceptionIsThrown() {
      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> configuration(List.of(), 100));
      assertThrows(IllegalArgumentException.class, () -> configuration(List.of("mafft"), 0));
    }

    private RunConfiguration configuration(List<String> seeds, int progressInterval) {
      return new RunConfiguration(
          directory,
          "classF",
          seeds,
          directory.resolve("out"),
          1L,
          AlgorithmConfiguration.withDefaults(100, 10, 1),
          MappingPolicy.VERIFIED,
          15f,
          progressInterval);
    }
  }

  @Nested
  @DisplayName("When locating precomputed alignments")
  class Precomputed {

    @ParameterizedTest
    @CsvSource({
      "mafft_classF.fasta, classF, mafft",
      "usalign_classF_clean.fasta, classF, usalign_clean",
      "MUSCLE5_CLASSF.FASTA, classF, MUSCLE5",
      "other.fasta, classF, other"
    })
    @DisplayName("given a file name, when deriving the method, then the dataset token is removed")
    void givenFileName_whenDerivingMethod_thenDatasetTokenIsRemoved(
        String file, String dataset, String method) {
      // Act & Assert
      assertEquals(method, PrecomputedAlignments.methodName(Path.of(file), dataset));
    }

    @Test
    @DisplayName("given a dataset directory, when finding methods, then unique matches are found")
    void givenDatasetDirectory_whenFindingMethods_thenUniqueMatchesAreFound() throws IOException {
      // Arrange
      Path precomputed = SyntheticDataset.write(directory).resolve("precomputed/synthetic");

      // Act & Assert
      assertEquals(
          Optional.of(precomputed.resolve("left_synthetic.fasta")),
          PrecomputedAlignments.find(precomputed, "synthetic", "left"));
      assertEquals(
          Optional.empty(), PrecomputedAlignments.find(precomputed, "synthetic", "missing"));
      assertEquals(
          Optional.empty(), PrecomputedAlignments.find(directory.resolve("none"), "x", "left"));
      assertEquals(3, PrecomputedAlignments.list(precomputed).size());
    }
  }

  @Nested
  @DisplayName("When recording progress")
  class Progress {

    @Test
    @DisplayName("given progress reports, when the interval is reached, then a line is appended")
    void givenProgressReports_whenIntervalIsReached_thenLineIsAppended() throws IOException {
      // Arrange
      ProteinSet proteins = Alignments.proteins("ACD:OOO", "AD:OO");
      StructuralMsaProblem problem =
          new StructuralMsaProblem(
              "p", proteins, List.of(TopologyAwareSumOfPairs.withDefaultParameters()));
      MsaSolution solution = Alignments.solution(proteins, "ACD", "A-D");
      problem.evaluate(solution);
      ProgressObserver observer = new ProgressObserver(problem, directory, 10);
      List<MsaSolution> front = new ArrayList<>(List.of(solution));

      // Act: 5 is below the interval, 15 follows the record at 12, the empty front is skipped
      observer.progress(5, () -> front);
      observer.progress(12, () -> front);
      observer.progress(15, () -> front);
      observer.progress(20, List::of);
      observer.progress(30, () -> front);

      // Assert
      List<String> lines = Files.readAllLines(directory.resolve(ProgressObserver.FILE_NAME));
      assertEquals(2, lines.size());
      assertTrue(lines.get(0).startsWith("{\"evaluations\":12,"));
      assertTrue(lines.get(1).contains("\"best\":{\"SOPwTP\":"));
    }

    @Test
    @DisplayName("given build information, when describing, then the expected keys are present")
    void givenBuildInformation_whenDescribing_thenExpectedKeysArePresent() {
      // Act
      Map<String, Object> info = BuildInfo.describe();

      // Assert
      assertEquals(List.of("gitCommit", "gitBranch", "gitDirty", "javaVersion"),
          List.copyOf(info.keySet()));
    }
  }
}
