package org.transept.seed;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.transept.alignment.GapColumnCompactor;
import org.transept.alignment.GapIntervals;
import org.transept.model.ProteinSet;
import org.transept.testing.Alignments;

@DisplayName("Unit tests for class InternalSeeds")
class InternalSeedsTest {

  private ProteinSet proteins;

  @BeforeEach
  void setUp() {
    proteins =
        Alignments.proteins(
            "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ:OOOOMMMMMMMMMMOOOOOOMMMMMMMMMOOOO",
            "MKTAYIAKQRISFVKSHFSRQLERLGLIEVQ:OOOOMMMMMMMMMOOOOOOMMMMMMMMOOOO",
            "MKAYIAKQRQISFVKHFSRQLEERLGIEVQA:OOOMMMMMMMMMMOOOOOMMMMMMMMMOOOO",
            "MTAYIAKQRQIFVKSHFRQLEERLGLIEQ:OOOMMMMMMMMMOOOOOMMMMMMMMMOOO");
  }

  @ParameterizedTest
  @ValueSource(ints = {1, 3, 6})
  @DisplayName("given a number of seeds, when generating, then the seeds are distinct and valid")
  void givenNumberOfSeeds_whenGenerating_thenTheSeedsAreDistinctAndValid(int count) {
    // Act
    List<InternalSeeds.Seed> seeds = InternalSeeds.generate(proteins, count);

    // Assert: at most the requested number, all different, every row keeps its residues
    assertFalse(seeds.isEmpty());
    assertTrue(seeds.size() <= count);
    Set<List<GapIntervals>> distinct = new HashSet<>();
    for (InternalSeeds.Seed seed : seeds) {
      assertTrue(distinct.add(GapColumnCompactor.compact(seed.alignment().rows())));
      assertTrue(seed.method().startsWith("biojava-g"));
      for (int row = 0; row < proteins.size(); row++) {
        assertEquals(
            proteins.get(row).residues(), seed.alignment().alignedRow(row).replace("-", ""));
      }
    }
    assertEquals(InternalSeeds.SETTINGS.get(0).method(), seeds.get(0).method());
  }

  @Test
  @DisplayName("given the same proteins, when generating twice, then the seeds are the same")
  void givenSameProteins_whenGeneratingTwice_thenTheSeedsAreTheSame() {
    // Act
    List<InternalSeeds.Seed> first = InternalSeeds.generate(proteins, 4);
    List<InternalSeeds.Seed> second = InternalSeeds.generate(proteins, 4);

    // Assert
    assertEquals(
        first.stream().map(seed -> seed.alignment().rows()).toList(),
        second.stream().map(seed -> seed.alignment().rows()).toList());
  }

  @Test
  @DisplayName("given an invalid number of seeds, when generating, then an exception is thrown")
  void givenInvalidNumberOfSeeds_whenGenerating_thenIllegalArgumentExceptionIsThrown() {
    // Act & Assert
    assertThrows(IllegalArgumentException.class, () -> InternalSeeds.generate(proteins, 0));
    assertThrows(
        IllegalArgumentException.class,
        () -> InternalSeeds.generate(proteins, InternalSeeds.MAX_SEEDS + 1));
  }
}
