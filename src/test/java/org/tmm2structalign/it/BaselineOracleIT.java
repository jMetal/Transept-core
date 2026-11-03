package org.tmm2structalign.it;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assumptions.assumeTrue;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.tmm2structalign.objective.LddtAggregation;
import org.tmm2structalign.run.BaselineScorer;
import org.tmm2structalign.run.BaselineScorer.BaselineScore;
import org.tmm2structalign.structure.MappingPolicy;

/**
 * Differential oracle: the objectives of the new core must reproduce the values that the
 * original implementation wrote to {@code precomputed/baselines.csv} (dataset v1.0.0, frozen in
 * {@code src/test/resources/oracle}).
 *
 * <p>Needs the real dataset; set {@code GPCR_DATASET} to the {@code GPCRdb} directory.
 */
@DisplayName("Integration tests: differential oracle against the original baselines.csv")
class BaselineOracleIT {

  /** Datasets whose lDDT changes under the VERIFIED policy (they contain OXER1 or MCHR1). */
  private static final Set<String> REMAPPED_DATASETS =
      Set.of("classA_002", "classA_004", "classA_004_19");

  /** Precision of the oracle file, which stores 10 decimals. */
  private static final double LDDT_TOLERANCE = 1e-9;

  private Path dataRoot;
  private List<String[]> expectedRows;

  @BeforeEach
  void setUp() throws IOException {
    String dataset = System.getenv("GPCR_DATASET");
    assumeTrue(dataset != null && Files.isDirectory(Path.of(dataset)), "GPCR_DATASET not set");
    dataRoot = Path.of(dataset);
    expectedRows = readOracle();
  }

  @Test
  @DisplayName("given the legacy mapping, when scoring baselines, then every value matches")
  void givenLegacyMapping_whenScoringBaselines_thenEveryValueMatches() throws IOException {
    // Act
    Map<String, BaselineScore> scores = score(MappingPolicy.LEGACY);

    // Assert
    int compared = 0;
    for (String[] row : expectedRows) {
      BaselineScore score = scores.get(row[0] + "/" + row[1]);
      assertEquals(Double.parseDouble(row[2]), score.values().get(0), row[0] + " " + row[1]);
      assertEquals(
          Double.parseDouble(row[3]), score.values().get(1), LDDT_TOLERANCE, row[0] + " " + row[1]);
      compared++;
    }
    assertEquals(162, compared);
  }

  @Test
  @DisplayName("given the verified mapping, when scoring baselines, then only remapped differ")
  void givenVerifiedMapping_whenScoringBaselines_thenOnlyRemappedDatasetsDiffer()
      throws IOException {
    // Act
    Map<String, BaselineScore> scores = score(MappingPolicy.VERIFIED);

    // Assert
    for (String[] row : expectedRows) {
      BaselineScore score = scores.get(row[0] + "/" + row[1]);
      String label = row[0] + " " + row[1];
      assertEquals(Double.parseDouble(row[2]), score.values().get(0), label);
      double difference = Math.abs(Double.parseDouble(row[3]) - score.values().get(1));
      if (REMAPPED_DATASETS.contains(row[0])) {
        assertTrue(difference > LDDT_TOLERANCE && difference < 0.02, label + ": " + difference);
      } else {
        assertTrue(difference <= LDDT_TOLERANCE, label + ": " + difference);
      }
    }
  }

  private Map<String, BaselineScore> score(MappingPolicy policy) throws IOException {
    BaselineScorer scorer =
        new BaselineScorer(dataRoot, policy, LddtAggregation.HARMONIC_OVER_LENGTH, 15f);
    Map<String, BaselineScore> scores = new HashMap<>();
    for (String dataset : expectedRows.stream().map(row -> row[0]).distinct().toList()) {
      for (BaselineScore score : scorer.score(dataset)) {
        scores.put(dataset + "/" + score.method(), score);
      }
    }
    return scores;
  }

  /** Returns dataset, method, objective 1 and objective 2 of every row with status "ok". */
  private static List<String[]> readOracle() throws IOException {
    List<String[]> rows = new ArrayList<>();
    try (InputStream input =
        BaselineOracleIT.class.getResourceAsStream("/oracle/baselines-v1.0.0.csv")) {
      String[] lines = new String(input.readAllBytes(), StandardCharsets.UTF_8).split("\n");
      for (int line = 1; line < lines.length; line++) {
        String[] fields = lines[line].split(",(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)", -1);
        if ("ok".equals(fields[5].trim())) {
          rows.add(
              new String[] {
                fields[0].replace("\"", ""), fields[1].replace("\"", ""), fields[2], fields[3]
              });
        }
      }
    }
    return rows;
  }
}
