package org.transept.run;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Instant;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import org.transept.algorithm.AlgorithmFactory;
import org.transept.algorithm.SearchRun;
import org.transept.algorithm.StopRequest;
import org.transept.alignment.GapIntervals;
import org.transept.alignment.MsaSolution;
import org.transept.io.FastaAlignmentWriter;
import org.transept.objective.Objective;
import org.transept.operator.RandomSource;
import org.transept.problem.StructuralMsaProblem;
import org.transept.structure.StructureMapping;
import org.uma.jmetal.util.pseudorandom.JMetalRandom;

/**
 * Runs one experiment and writes its results.
 *
 * <p>Output files in the run directory:
 *
 * <ul>
 *   <li>{@code FUN.csv}: one line per front solution, objective values in their natural
 *       direction, comma-separated, no header (jMetal convention; objective names are in
 *       {@code metadata.json})
 *   <li>{@code MSASol<i>.fasta}: alignment of the solution in line {@code i} of {@code FUN.csv}
 *   <li>{@code runtime.csv}: header {@code runtime_ms,runtime_seconds,evaluations} and one row
 *   <li>{@code seed.txt}: the random seed
 *   <li>{@code config.json}: the run configuration
 *   <li>{@code metadata.json}: build, timing, objectives, seed provenance and structure mapping
 *   <li>{@code progress.jsonl}: progress records written during the run
 * </ul>
 */
public final class ExperimentRunner {

  /** Name of the tool in metadata.json. */
  public static final String TOOL = "Transept";

  /** Objective values of the final solutions. */
  public static final String FUN_FILE = "FUN.csv";

  /** Runtime and number of evaluations. */
  public static final String RUNTIME_FILE = "runtime.csv";

  /**
   * Summary of a finished run.
   *
   * @param runDirectory directory with the results
   * @param evaluations evaluations performed
   * @param runtimeMillis wall-clock time of the search
   * @param frontSize number of solutions written
   */
  public record RunSummary(Path runDirectory, int evaluations, long runtimeMillis, int frontSize) {}

  /**
   * Runs an experiment.
   *
   * @param configuration run configuration
   * @return the summary of the run
   * @throws IOException if inputs cannot be read or outputs cannot be written
   * @throws RunDirectoryExistsException if the run directory already exists
   */
  public RunSummary run(RunConfiguration configuration) throws IOException {
    Path directory = configuration.runDirectory();
    if (Files.exists(directory)) {
      throw new RunDirectoryExistsException(directory);
    }
    Files.createDirectories(directory);
    write(directory.resolve("config.json"), Json.write(configuration.toMap()) + "\n");
    write(directory.resolve("seed.txt"), configuration.randomSeed() + System.lineSeparator());
    final Instant start = Instant.now();

    JMetalRandom.getInstance().setSeed(configuration.randomSeed());
    Experiment experiment = Experiment.load(configuration);
    StructuralMsaProblem problem = experiment.problem();
    StopRequest stop = StopRequest.byFile(directory.resolve(StopRequest.FILE_NAME));
    SearchRun search =
        new AlgorithmFactory(
                problem,
                experiment.seedAlignments(),
                configuration.algorithm(),
                experiment.structures(),
                RandomSource.jmetal())
            .stoppableBy(stop)
            .create(configuration.algorithmType());
    ProgressObserver progress =
        new ProgressObserver(problem, directory, configuration.progressInterval());
    search.onProgress(evaluations -> progress.progress(evaluations, search::front));
    search.run();

    List<MsaSolution> front = search.front();
    writeFront(problem, front, directory);
    long runtime = search.runtimeMillis();
    writeRuntime(directory, runtime, search.evaluations());
    RunSummary summary = new RunSummary(directory, search.evaluations(), runtime, front.size());
    Map<String, Object> metadata = metadata(experiment, configuration, search, summary, start);
    // Stopped by the user: the budget was not spent, and the run cannot be repeated exactly
    metadata.put("stoppedByUser", stop.wasRequested());
    write(directory.resolve("metadata.json"), Json.write(metadata) + "\n");
    Files.deleteIfExists(directory.resolve(StopRequest.FILE_NAME));
    return summary;
  }

  private static void writeFront(StructuralMsaProblem problem, List<MsaSolution> front, Path dir)
      throws IOException {
    StringBuilder fun = new StringBuilder();
    FastaAlignmentWriter writer = new FastaAlignmentWriter();
    for (int index = 0; index < front.size(); index++) {
      MsaSolution solution = front.get(index);
      List<String> values = new ArrayList<>();
      for (int objective = 0; objective < problem.numberOfObjectives(); objective++) {
        values.add(Double.toString(problem.naturalValue(solution, objective)));
      }
      fun.append(String.join(",", values)).append('\n');
      writer.write(solution.toMsa(), dir.resolve("MSASol" + index + ".fasta"));
    }
    write(dir.resolve(FUN_FILE), fun.toString());
  }

  private static void writeRuntime(Path directory, long runtimeMillis, int evaluations)
      throws IOException {
    write(
        directory.resolve(RUNTIME_FILE),
        "runtime_ms,runtime_seconds,evaluations\n"
            + runtimeMillis + ","
            + String.format(Locale.ROOT, "%.3f", runtimeMillis / 1000.0) + ","
            + evaluations + "\n");
  }

  private static Map<String, Object> metadata(
      Experiment experiment,
      RunConfiguration configuration,
      SearchRun search,
      RunSummary summary,
      Instant start) {
    Map<String, Object> metadata = new LinkedHashMap<>();
    metadata.put("tool", TOOL);
    metadata.put("algorithmType", configuration.algorithmType());
    metadata.put("algorithmName", search.name());
    metadata.put("configuration", search.configuration());
    metadata.put("build", BuildInfo.describe());
    metadata.put("startedAt", start.toString());
    metadata.put("finishedAt", Instant.now().toString());
    metadata.put("runtimeMs", summary.runtimeMillis());
    metadata.put("evaluations", summary.evaluations());
    metadata.put("frontSize", summary.frontSize());
    metadata.put(
        "objectives",
        experiment.problem().objectives().stream().map(Objective::name).toList());
    metadata.put("seeds", seedProvenance(experiment));
    metadata.put("structures", structureReport(experiment));
    metadata.put("proteinsWithoutStructure", experiment.structures().proteinsWithoutStructure());
    return metadata;
  }

  private static List<Map<String, Object>> seedProvenance(Experiment experiment) {
    Map<List<GapIntervals>, String> firstMethodOfAlignment = new HashMap<>();
    List<Map<String, Object>> seeds = new ArrayList<>();
    for (Experiment.SeedAlignment seed : experiment.seeds()) {
      String original = firstMethodOfAlignment.putIfAbsent(seed.alignment().rows(), seed.method());
      Map<String, Object> entry = new LinkedHashMap<>();
      entry.put("method", seed.method());
      entry.put("file", seed.file().toAbsolutePath().toString());
      entry.put("sha256", seed.sha256());
      entry.put("duplicateOf", original);
      entry.put("repairedGapRuns", seed.repairedGapRuns());
      entry.put("origin", seed.origin());
      seeds.add(entry);
    }
    return seeds;
  }

  private static List<Map<String, Object>> structureReport(Experiment experiment) {
    List<Map<String, Object>> report = new ArrayList<>();
    for (StructureMapping mapping : experiment.structures().mappings()) {
      Map<String, Object> entry = new LinkedHashMap<>();
      entry.put("protein", mapping.proteinId());
      entry.put("status", mapping.status());
      entry.put("identity", mapping.identity());
      entry.put("coverage", mapping.coverage());
      entry.put("detail", mapping.detail());
      report.add(entry);
    }
    return report;
  }

  private static void write(Path file, String content) throws IOException {
    Files.writeString(file, content, StandardCharsets.UTF_8);
  }
}
