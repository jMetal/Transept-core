package org.transept.algorithm;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicLong;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

@DisplayName("Unit tests for classes StopRequest and StoppableTermination")
class StopRequestTest {

  @TempDir Path directory;

  @Test
  @DisplayName("given the file appears, when checking, then it is seen only after the interval")
  void givenTheFileAppears_whenChecking_thenItIsSeenOnlyAfterTheInterval() throws IOException {
    // Arrange: a clock that only moves when the test says so
    AtomicLong clock = new AtomicLong();
    Path file = directory.resolve(StopRequest.FILE_NAME);
    StopRequest stop = new StopRequest(file, clock::get);
    assertFalse(stop.isRequested());

    // Act: the file appears, but the last look was just now
    Files.writeString(file, "");
    boolean tooSoon = stop.isRequested();
    clock.addAndGet(StopRequest.CHECK_INTERVAL_MILLIS * 1_000_000L);

    // Assert
    assertFalse(tooSoon);
    assertTrue(stop.isRequested());
    Files.delete(file);
    assertTrue(stop.isRequested(), "once requested, it stays requested");
    assertTrue(stop.wasRequested());
  }

  @Test
  @DisplayName("given no file, when checking, then the stop is never requested")
  void givenNoFile_whenChecking_thenTheStopIsNeverRequested() {
    assertFalse(StopRequest.NEVER.isRequested());
    assertFalse(StopRequest.byFile(directory.resolve("absent")).isRequested());
  }

  @Test
  @DisplayName("given a stop request, when testing the termination, then either condition ends it")
  void givenStopRequest_whenTestingTheTermination_thenEitherConditionEndsIt() throws IOException {
    // Arrange
    Path file = directory.resolve(StopRequest.FILE_NAME);
    StoppableTermination budgetLeft =
        new StoppableTermination(status -> false, StopRequest.byFile(file));
    StoppableTermination budgetSpent =
        new StoppableTermination(status -> true, StopRequest.NEVER);

    // Act & Assert
    assertFalse(budgetLeft.isMet(Map.of()));
    assertTrue(budgetSpent.isMet(Map.of()));
    Files.writeString(file, "");
    assertTrue(new StoppableTermination(status -> false, StopRequest.byFile(file)).isMet(Map.of()));
    assertThrows(IllegalArgumentException.class, () -> new StoppableTermination(null, null));
  }

  @Test
  @DisplayName("given an algorithm without results yet, when stopping, then it waits for them")
  void givenAlgorithmWithoutResultsYet_whenStopping_thenItWaitsForThem() {
    List<Integer> population = List.of(1, 2, 3, 4);
    assertFalse(
        StoppableTermination.hasResults(Map.of("POPULATION", population, "EVALUATIONS", 4)));
    assertTrue(
        StoppableTermination.hasResults(Map.of("POPULATION", population, "EVALUATIONS", 8)));
    assertFalse(
        StoppableTermination.hasResults(Map.of("POPULATION", List.of(), "EVALUATIONS", 0)));
    assertFalse(StoppableTermination.hasResults(Map.of("BEST_SOLUTIONS", List.of())));
    assertTrue(StoppableTermination.hasResults(Map.of("BEST_SOLUTIONS", List.of(1))));
    assertTrue(StoppableTermination.hasResults(Map.of()));
  }
}
