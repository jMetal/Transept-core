package org.tmm2structalign.run;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;
import org.tmm2structalign.alignment.MsaSolution;
import org.tmm2structalign.problem.StructuralMsaProblem;
import org.uma.jmetal.util.observable.Observable;
import org.uma.jmetal.util.observer.Observer;

/**
 * Records the progress of a run in {@code progress.jsonl}: one JSON object per line with the
 * number of evaluations, the elapsed time and the best value of every objective in the current
 * population. Designed to be followed by a GUI while the run is active.
 */
public final class ProgressObserver implements Observer<Map<String, Object>> {

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

  @Override
  @SuppressWarnings("unchecked")
  public void update(Observable<Map<String, Object>> observable, Map<String, Object> data) {
    Object evaluations = data.get("EVALUATIONS");
    Object population = data.get("POPULATION");
    if (evaluations instanceof Integer count
        && population instanceof List<?> solutions
        && count >= nextRecord) {
      record(count, (List<MsaSolution>) solutions);
      nextRecord = (count / interval + 1) * interval;
    }
  }

  /**
   * Appends one progress record.
   *
   * @param evaluations evaluations performed so far
   * @param population current population, evaluated
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
