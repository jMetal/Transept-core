package org.transept.evolver;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Map;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

@DisplayName("Unit tests for class MsaConfigurations")
class MsaConfigurationsTest {

  @Nested
  @DisplayName("When parsing and formatting configurations")
  class Parsing {

    @Test
    @DisplayName("given a configuration, when parsed and formatted, then it is unchanged")
    void givenConfiguration_whenParsedAndFormatted_thenItIsUnchanged() {
      // Arrange
      String text = "--selection tournament --selectionTournamentSize 3";

      // Act
      Map<String, String> values = MsaConfigurations.parse("  " + text + "\n");

      // Assert
      assertEquals(Map.of("selection", "tournament", "selectionTournamentSize", "3"), values);
      assertEquals(text, MsaConfigurations.format(values));
      assertArrayEquals(text.split(" "), MsaConfigurations.arguments(values));
      assertTrue(MsaConfigurations.parse(" ").isEmpty());
    }

    @ParameterizedTest
    @ValueSource(strings = {"--selection", "selection tournament", "--a 1 --a 2", "--a --b"})
    @DisplayName("given a malformed configuration, when parsed, then an exception is thrown")
    void givenMalformedConfiguration_whenParsed_thenIllegalArgumentExceptionIsThrown(String text) {
      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> MsaConfigurations.parse(text));
    }
  }

  @Nested
  @DisplayName("When reading the parameter space of NSGA-II")
  class NsgaTwo {

    @Test
    @DisplayName("given the default configuration, when read, then it is the published one")
    void givenDefaultConfiguration_whenRead_thenItIsThePublishedOne() {
      // Act
      Map<String, String> defaults = MsaConfigurations.defaults(MsaNsga2.PARAMETER_SPACE);

      // Assert
      assertEquals("0.8", defaults.get("crossoverProbability"));
      assertEquals("0.2", defaults.get("mutationProbability"));
      assertEquals(MsaOperators.TM_SAFE_GAP_SHIFT, defaults.get("mutation"));
      assertEquals("2", defaults.get("selectionTournamentSize"));
      assertEquals("0", defaults.get("maxShift"));
      assertEquals("1", defaults.get("cutCandidates"));
      MsaConfigurations.requireKnownNames(
          defaults, MsaConfigurations.parameterSpace(MsaNsga2.PARAMETER_SPACE));
    }

    @Test
    @DisplayName("given an unknown parameter name, when checked, then an exception is thrown")
    void givenUnknownParameterName_whenChecked_thenIllegalArgumentExceptionIsThrown() {
      // Arrange
      var space = MsaConfigurations.parameterSpace(MsaNsga2.PARAMETER_SPACE);

      // Act & Assert
      assertThrows(
          IllegalArgumentException.class,
          () -> MsaConfigurations.requireKnownNames(Map.of("mutationProbabilty", "0.1"), space));
    }

    @Test
    @DisplayName("given a space without default configuration, when read, then it fails")
    void givenSpaceWithoutDefault_whenRead_thenIllegalArgumentExceptionIsThrown() {
      // Act & Assert
      assertThrows(
          IllegalArgumentException.class, () -> MsaConfigurations.defaults("Missing.yaml"));
    }
  }
}
