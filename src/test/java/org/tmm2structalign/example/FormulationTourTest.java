package org.tmm2structalign.example;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

@DisplayName("Unit tests for class FormulationTour")
class FormulationTourTest {

  private static final Path TUTORIAL = Path.of("docs/tutorials/T2-problem-formulation.md");

  @Test
  @DisplayName("given the tutorial, when running the tour, then every section is shown unchanged")
  void givenTutorial_whenRunningTour_thenEverySectionIsShownUnchanged() throws IOException {
    // Arrange
    String tutorial = Files.readString(TUTORIAL);

    // Act
    List<String> sections = sections(FormulationTour.run());

    // Assert
    assertEquals(6, sections.size());
    for (String section : sections) {
      assertTrue(
          tutorial.contains(section),
          "T2 does not show the current output of the tour; update it with:\n" + section);
    }
  }

  @Test
  @DisplayName("given the same seed, when running twice, then the output is identical")
  void givenSameSeed_whenRunningTwice_thenOutputIsIdentical() {
    assertEquals(FormulationTour.run(), FormulationTour.run());
  }

  @Test
  @DisplayName("given the main method, when called, then it prints the tour")
  void givenMainMethod_whenCalled_thenItPrintsTheTour() {
    // Arrange
    PrintStream original = System.out;
    ByteArrayOutputStream buffer = new ByteArrayOutputStream();

    // Act
    try {
      System.setOut(new PrintStream(buffer, true, StandardCharsets.UTF_8));
      FormulationTour.main(new String[0]);
    } finally {
      System.setOut(original);
    }

    // Assert
    assertEquals(FormulationTour.run(), buffer.toString(StandardCharsets.UTF_8));
  }

  /** Splits the output at the "== title ==" lines, keeping each heading with its section. */
  private static List<String> sections(String output) {
    List<String> sections = new ArrayList<>();
    StringBuilder current = new StringBuilder();
    for (String line : output.split("\n")) {
      if (line.startsWith("== ") && current.length() > 0) {
        sections.add(current.toString());
        current.setLength(0);
      }
      current.append(line).append('\n');
    }
    sections.add(current.toString());
    return sections;
  }
}
