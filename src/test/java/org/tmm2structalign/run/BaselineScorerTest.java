package org.tmm2structalign.run;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.tmm2structalign.io.InvalidInputFileException;
import org.tmm2structalign.objective.LddtAggregation;
import org.tmm2structalign.run.BaselineScorer.BaselineScore;
import org.tmm2structalign.structure.MappingPolicy;
import org.tmm2structalign.testing.SyntheticDataset;

@DisplayName("Unit tests for class BaselineScorer")
class BaselineScorerTest {

  @TempDir Path root;

  private BaselineScorer scorer;

  @BeforeEach
  void setUp() throws IOException {
    SyntheticDataset.write(root);
    Files.writeString(
        root.resolve("precomputed/synthetic/broken_synthetic.fasta"), ">prot0_synth\nAAAA\n");
    scorer =
        new BaselineScorer(root, MappingPolicy.VERIFIED, LddtAggregation.HARMONIC_OVER_LENGTH, 15f);
  }

  @Nested
  @DisplayName("When scoring a dataset")
  class Scoring {

    @Test
    @DisplayName("given seed files, when scoring, then each file is scored or reported")
    void givenSeedFiles_whenScoring_thenEachFileIsScoredOrReported() throws IOException {
      // Act
      List<BaselineScore> scores = scorer.score(SyntheticDataset.NAME);

      // Assert
      assertEquals(4, scores.size());
      BaselineScore broken = scores.get(0);
      assertEquals("broken", broken.method());
      assertFalse(broken.isScored());
      assertTrue(broken.values().isEmpty());
      BaselineScore trueAlignment = scores.get(3);
      assertEquals("true", trueAlignment.method());
      assertEquals("precomputed/synthetic/true_synthetic.fasta", trueAlignment.file());
      assertTrue(trueAlignment.isScored());
      assertEquals(2, trueAlignment.values().size());
      assertEquals(List.of("SOPwTP", "lDDT"), BaselineScorer.objectiveNames());
    }

    @Test
    @DisplayName("given the true alignment, when scoring, then it beats the padded ones")
    void givenTrueAlignment_whenScoring_thenItBeatsPaddedOnes() throws IOException {
      // Act
      List<BaselineScore> scores = scorer.score(SyntheticDataset.NAME);

      // Assert
      double trueSop = scores.get(3).values().get(0);
      double trueLddt = scores.get(3).values().get(1);
      for (BaselineScore padded : List.of(scores.get(1), scores.get(2))) {
        assertTrue(trueSop > padded.values().get(0), padded.method());
        assertTrue(trueLddt > padded.values().get(1), padded.method());
      }
    }

    @Test
    @DisplayName("given a missing dataset, when scoring, then an exception is thrown")
    void givenMissingDataset_whenScoring_thenInvalidInputFileExceptionIsThrown() {
      // Act & Assert
      assertThrows(InvalidInputFileException.class, () -> scorer.score("missing"));
    }
  }
}
