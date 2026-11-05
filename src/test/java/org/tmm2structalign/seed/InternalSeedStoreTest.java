package org.tmm2structalign.seed;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.tmm2structalign.model.ProteinSet;
import org.tmm2structalign.testing.Alignments;

@DisplayName("Unit tests for class InternalSeedStore")
class InternalSeedStoreTest {

  @TempDir Path directory;

  private ProteinSet proteins;
  private List<InternalSeeds.Seed> seeds;

  @BeforeEach
  void setUp() {
    proteins =
        Alignments.proteins(
            "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ:OOOOMMMMMMMMMMOOOOOOMMMMMMMMMOOOO",
            "MKTAYIAKQRISFVKSHFSRQLERLGLIEVQ:OOOOMMMMMMMMMOOOOOOMMMMMMMMOOOO",
            "MKAYIAKQRQISFVKHFSRQLEERLGIEVQA:OOOMMMMMMMMMMOOOOOMMMMMMMMMOOOO");
    seeds = InternalSeeds.generate(proteins, 3);
  }

  @Test
  @DisplayName("given written seeds, when read with the same proteins, then they are returned")
  void givenWrittenSeeds_whenReadWithSameProteins_thenTheyAreReturned() throws IOException {
    // Arrange
    InternalSeedStore.write(directory, "toy", proteins, seeds, "test");

    // Act
    InternalSeedStore.Lookup lookup = InternalSeedStore.read(directory, proteins);

    // Assert
    assertTrue(lookup.usable(), lookup.problem());
    assertEquals(
        seeds.stream().map(InternalSeeds.Seed::method).toList(),
        lookup.seeds().stream().map(InternalSeedStore.StoredSeed::method).toList());
    assertEquals(seeds.get(0).alignment().rows(), lookup.seeds().get(0).alignment().rows());
  }

  @Test
  @DisplayName("given other proteins, when reading, then the seeds are not used")
  void givenOtherProteins_whenReading_thenTheSeedsAreNotUsed() throws IOException {
    // Arrange: the same residues with another topology
    InternalSeedStore.write(directory, "toy", proteins, seeds, "test");
    ProteinSet other =
        Alignments.proteins(
            "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ:OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
            "MKTAYIAKQRISFVKSHFSRQLERLGLIEVQ:OOOOMMMMMMMMMOOOOOOMMMMMMMMOOOO",
            "MKAYIAKQRQISFVKHFSRQLEERLGIEVQA:OOOMMMMMMMMMMOOOOOMMMMMMMMMOOOO");

    // Act
    InternalSeedStore.Lookup lookup = InternalSeedStore.read(directory, other);

    // Assert
    assertFalse(lookup.usable());
    assertTrue(lookup.problem().contains("other proteins"), lookup.problem());
  }

  @Test
  @DisplayName("given a changed or missing file, when reading, then the seeds are not used")
  void givenChangedOrMissingFile_whenReading_thenTheSeedsAreNotUsed() throws IOException {
    // Arrange
    InternalSeedStore.write(directory, "toy", proteins, seeds, "test");
    Path first = directory.resolve(seeds.get(0).method() + "_toy.fasta");
    Files.writeString(first, Files.readString(first) + "\n");

    // Act & Assert
    assertFalse(InternalSeedStore.read(directory, proteins).usable());
    assertFalse(InternalSeedStore.read(directory.resolve("none"), proteins).usable());
  }
}
