package org.transept.cli;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Duration;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import org.transept.algorithm.AlgorithmConfiguration;
import org.transept.algorithm.AlgorithmType;
import org.transept.algorithm.CrossoverType;
import org.transept.algorithm.MutationType;
import org.transept.evolver.MsaConfigurations;
import org.transept.objective.ObjectiveSet;
import org.transept.objective.StructuralLddt;
import org.transept.operator.TmSafeGapShiftMutation;
import org.transept.operator.TmSafeSinglePointCrossover;
import org.transept.run.ExperimentRunner;
import org.transept.run.ExperimentRunner.RunSummary;
import org.transept.run.Json;
import org.transept.run.RunConfiguration;
import org.transept.run.RunDirectoryExistsException;
import org.transept.seed.InternalSeedStore;
import org.transept.seed.InternalSeeds;
import org.transept.structure.MappingPolicy;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

/**
 * Runs the algorithm on one dataset. Results go to
 * {@code <output-root>/ejecuciones/<dataset>/<run-id>/}, the layout read by
 * {@code scripts/pareto_analysis.py}.
 */
@Command(
    name = "run",
    mixinStandardHelpOptions = true,
    description = "Align one dataset and write the final population.")
public final class RunCommand implements Callable<Integer> {

  /** Name of the error report written to the run directory when a run fails. */
  public static final String ERROR_FILE = "error.json";

  /** Default mutation probability of PAES: every iteration produces a new mutant. */
  public static final double PAES_MUTATION_PROBABILITY = 1.0;

  @Option(names = "--data-root", required = true, description = "Dataset root (.../GPCRdb).")
  Path dataRoot;

  @Option(names = "--dataset", required = true, description = "Dataset name, e.g. classF.")
  String dataset;

  @Option(
      names = "--seeds",
      split = ",",
      defaultValue = InternalSeeds.NAME,
      description =
          "Seed alignments of the population: tools whose precomputed alignments are in the"
              + " dataset, e.g. mafft,kalign, and/or internal (up to 10 alignments computed by the"
              + " tool with BioJava). Default: ${DEFAULT-VALUE}.")
  List<String> seeds;

  @Option(
      names = "--repair-seeds",
      negatable = true,
      defaultValue = "true",
      fallbackValue = "true",
      description =
          "Move the gaps that split a transmembrane helix out of the seed alignments before the"
              + " search (the operators never create such gaps). Default: yes."
              + " --no-repair-seeds keeps the gaps of the precomputed alignments; the internal"
              + " seeds are always repaired.")
  boolean repairSeeds;

  @Option(names = "--output-root", required = true, description = "Root of the results tree.")
  Path outputRoot;

  @Option(names = "--run-id", required = true, description = "Run identifier, e.g. seed_000.")
  String runId;

  @Option(names = "--seed", required = true, description = "Seed of the random generator.")
  long seed;

  /** Evaluations of a run when no stopping condition is given. */
  public static final int DEFAULT_MAX_EVALUATIONS = 10_000;

  @Option(
      names = "--max-evaluations",
      description = "Stop after this many evaluations. Default: 10000, unless --max-time is given.")
  Integer maxEvaluations;

  @Option(
      names = "--max-time",
      description =
          "Stop after this many seconds of search instead of a number of evaluations. Runs that"
              + " stop by time are not reproducible exactly.")
  Double maxTime;

  @Option(
      names = "--population-size",
      defaultValue = "100",
      description = "Default: ${DEFAULT-VALUE}.")
  int populationSize;

  @Option(
      names = "--threads",
      defaultValue = "1",
      description = "Evaluation threads. Default: ${DEFAULT-VALUE}.")
  int threads;

  @Option(names = "--crossover-probability", description = "Default: 0.8.")
  Double crossoverProbability;

  @Option(
      names = "--crossover",
      description =
          "Crossover operator: ${COMPLETION-CANDIDATES}. Default: TM_SAFE_SINGLE_POINT. Neither"
              + " cuts a transmembrane helix. TM_SAFE_UNIFORM exchanges blocks between the"
              + " columns where both parents agree, adding no gap.")
  CrossoverType crossover;

  @Option(
      names = "--cut-candidates",
      description =
          "Valid cuts the crossover draws, of which it uses the one that adds the fewest gaps."
              + " Default: 1, a cut at random.")
  Integer cutCandidates;

  @Option(
      names = "--mutation-probability",
      description = "Default: 1.0 for PAES (mutation-only), 0.2 for the other algorithms.")
  Double mutationProbability;

  @Option(
      names = "--mutation",
      description =
          "Mutation operator: ${COMPLETION-CANDIDATES}. Default: TM_SAFE_GAP_SHIFT. None of them"
              + " breaks a transmembrane helix. RANDOM_CHOICE picks one of the other four each"
              + " time.")
  MutationType mutation;

  @Option(
      names = "--max-shift",
      description =
          "Largest number of columns the mutation can move a gap run (TM_SAFE_GAP_SHIFT and"
              + " RANDOM_CHOICE). Default: 0, no limit.")
  Integer maxShift;

  @Option(
      names = "--algorithm",
      defaultValue = "NSGA_II",
      description =
          "Search algorithm: ${COMPLETION-CANDIDATES}. Default: ${DEFAULT-VALUE}. NSGA_III rounds"
              + " the population up to a multiple of 4; with three objectives, it and MOEA_D need a"
              + " population size with weight vectors; PAES uses the first seed and ignores the"
              + " crossover; RANDOM_SEARCH ignores the seeds and both operators. The ASYNC_"
              + " algorithms keep all --threads busy but do not repeat exactly.")
  AlgorithmType algorithm;

  @Option(
      names = "--configuration",
      description =
          "Configuration of NSGA_II, SPEA2 or GENETIC_ALGORITHM from its parameter space, as"
              + " '--parameter value ...' (see docs/reference/configuration.md). Unset"
              + " parameters keep their default value.")
  String configuration;

  @Option(
      names = "--configuration-file",
      description = "File whose first non-blank line is a configuration (see --configuration).")
  Path configurationFile;

  @Option(
      names = "--objectives",
      defaultValue = "SOPWTP_LDDT",
      description =
          "Objectives to optimise: ${COMPLETION-CANDIDATES}. Default: ${DEFAULT-VALUE}"
              + " (published). LDDT needs a single-objective algorithm (GENETIC_ALGORITHM,"
              + " ASYNC_GENETIC_ALGORITHM), and they need LDDT.")
  ObjectiveSet objectives;

  @Option(
      names = "--mapping-policy",
      defaultValue = "VERIFIED",
      description =
          "Sequence-structure mapping: ${COMPLETION-CANDIDATES}. Default: ${DEFAULT-VALUE}.")
  MappingPolicy mappingPolicy;

  @Option(
      names = "--progress-interval",
      defaultValue = "100",
      description = "Evaluations between progress records. Default: ${DEFAULT-VALUE}.")
  int progressInterval;

  @Override
  public Integer call() throws IOException {
    Path runDirectory = outputRoot.resolve("ejecuciones").resolve(dataset).resolve(runId);
    int status = 0;
    try {
      RunConfiguration configuration =
          new RunConfiguration(
              dataRoot,
              dataset,
              seeds,
              runDirectory,
              seed,
              new AlgorithmConfiguration(
                  maxEvaluations == null ? DEFAULT_MAX_EVALUATIONS : maxEvaluations,
                  populationSize,
                  populationSize,
                  crossoverProbability == null
                      ? AlgorithmConfiguration.DEFAULT_CROSSOVER_PROBABILITY
                      : crossoverProbability,
                  effectiveMutationProbability(),
                  threads,
                  mutation == null ? MutationType.TM_SAFE_GAP_SHIFT : mutation,
                  maxRunTime(),
                  configurationText(),
                  maxShift == null ? TmSafeGapShiftMutation.UNLIMITED_SHIFT : maxShift,
                  cutCandidates == null
                      ? TmSafeSinglePointCrossover.DEFAULT_CUT_CANDIDATES
                      : cutCandidates,
                  crossover == null ? CrossoverType.TM_SAFE_SINGLE_POINT : crossover),
              algorithm,
              mappingPolicy,
              StructuralLddt.DEFAULT_CUTOFF,
              progressInterval,
              objectives,
              repairSeeds,
              InternalSeedStore.directory(outputRoot, dataset));
      RunSummary summary = new ExperimentRunner().run(configuration);
      System.out.println(
          "Run finished: " + summary.evaluations() + " evaluations, " + summary.frontSize()
              + " solutions in " + summary.runDirectory());
    } catch (RunDirectoryExistsException exception) {
      System.err.println(exception.getMessage());
      status = 2;
    } catch (RuntimeException | IOException exception) {
      reportError(runDirectory, exception);
      status = 1;
    }
    return status;
  }

  /**
   * Configuration given with {@code --configuration} or {@code --configuration-file}, or an empty
   * string.
   *
   * @throws IllegalArgumentException if both are given, or if an operator option sets a parameter
   *     that the configuration also sets
   * @throws IOException if the file cannot be read
   */
  String configurationText() throws IOException {
    if (configuration != null && configurationFile != null) {
      throw new IllegalArgumentException(
          "Give --configuration or --configuration-file, not both");
    }
    String text = configuration == null ? "" : configuration;
    if (configurationFile != null) {
      text =
          Files.readAllLines(configurationFile, StandardCharsets.UTF_8).stream()
              .filter(line -> !line.isBlank())
              .findFirst()
              .orElse("");
    }
    Map<String, String> values = MsaConfigurations.parse(text);
    requireNotBoth(crossoverProbability != null, values, "crossoverProbability",
        "--crossover-probability");
    requireNotBoth(mutationProbability != null, values, "mutationProbability",
        "--mutation-probability");
    requireNotBoth(crossover != null, values, "crossover", "--crossover");
    requireNotBoth(mutation != null, values, "mutation", "--mutation");
    requireNotBoth(maxShift != null, values, "maxShift", "--max-shift");
    requireNotBoth(cutCandidates != null, values, "cutCandidates", "--cut-candidates");
    return text;
  }

  private static void requireNotBoth(
      boolean optionGiven, Map<String, String> values, String parameter, String option) {
    if (optionGiven && values.containsKey(parameter)) {
      throw new IllegalArgumentException(
          option + " and --" + parameter + " in the configuration set the same parameter");
    }
  }

  /**
   * Time limit given with {@code --max-time}, or {@code null} to stop by evaluations.
   *
   * @throws IllegalArgumentException if both stopping conditions are given or the time is not
   *     positive
   */
  Duration maxRunTime() {
    if (maxTime == null) {
      return null;
    }
    if (maxEvaluations != null) {
      throw new IllegalArgumentException("Give --max-evaluations or --max-time, not both");
    }
    if (!(maxTime > 0)) {
      throw new IllegalArgumentException("--max-time must be positive: " + maxTime);
    }
    return Duration.ofMillis(Math.round(maxTime * 1000));
  }

  /** Mutation probability given by the user, or the default of the chosen algorithm. */
  double effectiveMutationProbability() {
    double defaultValue =
        algorithm == AlgorithmType.PAES
            ? PAES_MUTATION_PROBABILITY
            : AlgorithmConfiguration.DEFAULT_MUTATION_PROBABILITY;
    return mutationProbability != null ? mutationProbability : defaultValue;
  }

  private static void reportError(Path runDirectory, Exception exception) throws IOException {
    Map<String, Object> error = new LinkedHashMap<>();
    error.put("type", exception.getClass().getName());
    error.put("message", String.valueOf(exception.getMessage()));
    String json = Json.write(error);
    System.err.println(json);
    if (Files.isDirectory(runDirectory)) {
      Files.writeString(
          runDirectory.resolve(ERROR_FILE), json + "\n", StandardCharsets.UTF_8);
    }
  }
}
