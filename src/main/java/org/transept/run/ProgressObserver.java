package org.transept.run;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Supplier;
import java.util.logging.Logger;
import org.transept.alignment.MsaSolution;
import org.transept.problem.StructuralMsaProblem;

/**
 * Records the progress of a run in {@code progress.jsonl}: one JSON object per line with the
 * number of evaluations, the elapsed time and the best value of every objective in the current
 * result of the search (population, archive or non-dominated solutions, depending on the
 * algorithm). Designed to be followed by a GUI while the run is active.
 */
public final class ProgressObserver {

  /** File name of the progress log inside the run directory. */
  public static final String FILE_NAME = "progress.jsonl";

  private static final Logger LOGGER = Logger.getLogger(ProgressObserver.class.getName());

  private final StructuralMsaProblem problem;
  private final Path file;
  private final int interval;
  private final long startNanos;
  private int nextRecord;

  /**
   * Creates the observer.
   *
   * @param problem problem being solved, used to name and orient objective values
   * @param runDirectory directory where the progress log is written
   * @param interval minimum number of evaluations between records
   */
  public ProgressObserver(StructuralMsaProblem problem, Path runDirectory, int interval) {
    this.problem = problem;
    this.file = runDirectory.resolve(FILE_NAME);
    this.interval = interval;
    this.startNanos = System.nanoTime();
    this.nextRecord = interval;
  }

  /**
   * Records progress if at least {@code interval} evaluations have passed since the last record.
   *
   * @param evaluations evaluations performed so far
   * @param front supplier of the current result of the search
   */
  public void progress(int evaluations, Supplier<List<MsaSolution>> front) {
    if (evaluations >= nextRecord) {
      List<MsaSolution> solutions = front.get();
      if (!solutions.isEmpty()) {
        record(evaluations, solutions);
      }
      nextRecord = (evaluations / interval + 1) * interval;
    }
  }

  /**
   * Appends one progress record.
   *
   * @param evaluations evaluations performed so far
   * @param population current result of the search, evaluated
   */
  public void record(int evaluations, List<MsaSolution> population) {
    Map<String, Object> best = new LinkedHashMap<>();
    for (int index = 0; index < problem.numberOfObjectives(); index++) {
      boolean maximized = problem.objectives().get(index).isMaximized();
      double bestValue = maximized ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
      for (MsaSolution solution : population) {
        double value = problem.naturalValue(solution, index);
        bestValue = maximized ? Math.max(bestValue, value) : Math.min(bestValue, value);
      }
      best.put(problem.objectives().get(index).name(), bestValue);
    }
    Map<String, Object> line = new LinkedHashMap<>();
    line.put("evaluations", evaluations);
    line.put("elapsedMs", (System.nanoTime() - startNanos) / 1_000_000L);
    line.put("best", best);
    String json = Json.write(line);
    LOGGER.fine(json);
    try {
      Files.writeString(
          file,
          json + "\n",
          StandardCharsets.UTF_8,
          StandardOpenOption.CREATE,
          StandardOpenOption.APPEND);
    } catch (IOException exception) {
      throw new UncheckedIOException("Cannot write progress to " + file, exception);
    }
  }
}
