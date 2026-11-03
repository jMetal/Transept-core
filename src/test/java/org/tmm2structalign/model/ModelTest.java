package org.tmm2structalign.model;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import java.util.Optional;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

@DisplayName("Unit tests for the model classes Topology, ProteinSequence and ProteinSet")
class ModelTest {

  @Nested
  @DisplayName("When decoding topology codes")
  class TopologyCodes {

    @ParameterizedTest
    @CsvSource({
      "M, TRANSMEMBRANE",
      "I, INSIDE",
      "O, OUTSIDE",
      "S, SIGNAL_PEPTIDE",
      "X, UNKNOWN",
      "m, UNKNOWN"
    })
    @DisplayName("given a code, when decoding, then the topology is returned")
    void givenCode_whenDecoding_thenTopologyIsReturned(char code, Topology expected) {
      // Act & Assert
      assertEquals(expected, Topology.fromCode(code));
    }

    @Test
    @DisplayName("given each topology, when reading its code, then it round-trips")
    void givenEachTopology_whenReadingCode_thenItRoundTrips() {
      // Act & Assert
      for (Topology topology : Topology.values()) {
        assertEquals(topology, Topology.fromCode(topology.code()));
      }
      assertTrue(Topology.TRANSMEMBRANE.isTransmembrane());
      assertFalse(Topology.OUTSIDE.isTransmembrane());
    }
  }

  @Nested
  @DisplayName("When creating protein sequences")
  class Sequences {

    @Test
    @DisplayName("given residues and topology, when reading, then positions are exposed")
    void givenResiduesAndTopology_whenReading_thenPositionsAreExposed() {
      // Act
      ProteinSequence protein = new ProteinSequence("p", "ACD", "OMI");

      // Assert
      assertEquals(3, protein.length());
      assertEquals('C', protein.residueAt(1));
      assertEquals(Topology.TRANSMEMBRANE, protein.topologyAt(1));
      assertTrue(protein.isTransmembrane(1));
      assertFalse(protein.isTransmembrane(2));
    }

    @Test
    @DisplayName("given no topology, when creating, then every residue is unknown")
    void givenNoTopology_whenCreating_thenEveryResidueIsUnknown() {
      // Act
      ProteinSequence protein = ProteinSequence.withoutTopology("p", "AC");

      // Assert
      assertEquals("UU", protein.topology());
    }

    @ParameterizedTest
    @CsvSource({"'', ACD, OOO", "' ', ACD, OOO", "p, '', ''", "p, A-C, OOO", "p, ACD, OO"})
    @DisplayName("given invalid fields, when creating, then an exception is thrown")
    void givenInvalidFields_whenCreating_thenIllegalArgumentExceptionIsThrown(
        String id, String residues, String topology) {
      // Act & Assert
      assertThrows(
          IllegalArgumentException.class, () -> new ProteinSequence(id, residues, topology));
    }
  }

  @Nested
  @DisplayName("When creating protein sets")
  class Sets {

    @Test
    @DisplayName("given distinct proteins, when looking up ids, then indices are returned")
    void givenDistinctProteins_whenLookingUpIds_thenIndicesAreReturned() {
      // Act
      ProteinSet set =
          new ProteinSet(
              List.of(
                  ProteinSequence.withoutTopology("a", "AC"),
                  ProteinSequence.withoutTopology("b", "DE")));

      // Assert
      assertEquals(2, set.size());
      assertEquals(Optional.of(1), set.indexOf("b"));
      assertEquals(Optional.empty(), set.indexOf("c"));
      assertEquals("a", set.get(0).id());
      assertEquals(2, set.proteins().size());
    }

    @Test
    @DisplayName("given duplicated ids, when creating, then an exception is thrown")
    void givenDuplicatedIds_whenCreating_thenIllegalArgumentExceptionIsThrown() {
      // Arrange
      List<ProteinSequence> proteins =
          List.of(
              ProteinSequence.withoutTopology("a", "AC"),
              ProteinSequence.withoutTopology("a", "DE"));

      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> new ProteinSet(proteins));
    }

    @Test
    @DisplayName("given fewer than two proteins, when creating, then an exception is thrown")
    void givenSingleProtein_whenCreating_thenIllegalArgumentExceptionIsThrown() {
      // Arrange
      List<ProteinSequence> proteins = List.of(ProteinSequence.withoutTopology("a", "AC"));

      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> new ProteinSet(proteins));
    }
  }
}
