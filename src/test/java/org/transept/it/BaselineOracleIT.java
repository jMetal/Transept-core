package org.transept.it;

import static org.junit.jupiter.api.Assertions.assertEquals;
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
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.transept.run.BaselineScorer;
import org.transept.run.BaselineScorer.BaselineScore;
import org.transept.structure.MappingPolicy;

/**
 * Differential oracle for the objectives on the precomputed alignments of dataset v1.0.0.
 *
 * <p>SOPwTP must reproduce the values that the original implementation wrote to {@code
 * precomputed/baselines.csv} (frozen in {@code src/test/resources/oracle}). The lDDT objective no
 * longer uses the original aggregation, so its values are compared with a regression snapshot
 * computed by this version ({@code lddt-snapshot.csv}); regenerate it only when the objective
 * changes on purpose, and explain why.
 *
 * <p>Needs the real dataset; set {@code GPCR_DATASET} to the {@code GPCRdb} directory.
 */
@DisplayName("Integration tests: differential oracle against the original baselines.csv")
class BaselineOracleIT {

  /** Precision of the snapshot, which stores 12 decimals. */
  private static final double LDDT_TOLERANCE = 1e-11;

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
  @DisplayName("given the baselines, when scoring, then SOPwTP matches the original everywhere")
  void givenBaselines_whenScoring_thenSopwtpMatchesTheOriginalEverywhere() throws IOException {
    // Act
    Map<String, BaselineScore> scores = score(MappingPolicy.VERIFIED);

    // Assert
    int compared = 0;
    for (String[] row : expectedRows) {
      BaselineScore score = scores.get(row[0] + "/" + row[1]);
      assertEquals(Double.parseDouble(row[2]), score.values().get(0), row[0] + " " + row[1]);
      compared++;
    }
    assertEquals(162, compared);
  }

  @Test
  @DisplayName("given the baselines, when scoring, then lDDT matches the regression snapshot")
  void givenBaselines_whenScoring_thenLddtMatchesTheRegressionSnapshot() throws IOException {
    // Arrange
    Map<String, Double> snapshot = readSnapshot();

    // Act
    Map<String, BaselineScore> scores = score(MappingPolicy.VERIFIED);

    // Assert
    assertEquals(162, snapshot.size());
    for (Map.Entry<String, Double> expected : snapshot.entrySet()) {
      BaselineScore score = scores.get(expected.getKey());
      assertEquals(expected.getValue(), score.values().get(1), LDDT_TOLERANCE, expected.getKey());
    }
  }

  private Map<String, BaselineScore> score(MappingPolicy policy) throws IOException {
    BaselineScorer scorer =
        new BaselineScorer(dataRoot, policy, 15f);
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

  /** Returns the snapshot lDDT of every precomputed alignment, keyed by dataset/method. */
  private static Map<String, Double> readSnapshot() throws IOException {
    Map<String, Double> snapshot = new HashMap<>();
    try (InputStream input =
        BaselineOracleIT.class.getResourceAsStream("/oracle/lddt-snapshot.csv")) {
      String[] lines = new String(input.readAllBytes(), StandardCharsets.UTF_8).split("\n");
      for (int line = 1; line < lines.length; line++) {
        String[] fields = lines[line].split(",");
        snapshot.put(fields[0] + "/" + fields[1], Double.parseDouble(fields[2]));
      }
    }
    return snapshot;
  }
}
