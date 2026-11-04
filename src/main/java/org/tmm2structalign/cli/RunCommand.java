package org.tmm2structalign.cli;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import org.tmm2structalign.algorithm.AlgorithmConfiguration;
import org.tmm2structalign.algorithm.AlgorithmType;
import org.tmm2structalign.algorithm.InitializationType;
import org.tmm2structalign.algorithm.MutationType;
import org.tmm2structalign.objective.StructuralLddt;
import org.tmm2structalign.run.ExperimentRunner;
import org.tmm2structalign.run.ExperimentRunner.RunSummary;
import org.tmm2structalign.run.Json;
import org.tmm2structalign.run.RunConfiguration;
import org.tmm2structalign.run.RunDirectoryExistsException;
import org.tmm2structalign.structure.MappingPolicy;
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
      required = true,
      split = ",",
      description = "Tools whose precomputed alignments seed the population, e.g. mafft,kalign.")
  List<String> seeds;

  @Option(names = "--output-root", required = true, description = "Root of the results tree.")
  Path outputRoot;

  @Option(names = "--run-id", required = true, description = "Run identifier, e.g. seed_000.")
  String runId;

  @Option(names = "--seed", required = true, description = "Seed of the random generator.")
  long seed;

  @Option(
      names = "--max-evaluations",
      defaultValue = "10000",
      description = "Default: ${DEFAULT-VALUE}.")
  int maxEvaluations;

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

  @Option(
      names = "--crossover-probability",
      defaultValue = "0.8",
      description = "Default: ${DEFAULT-VALUE}.")
  double crossoverProbability;

  @Option(
      names = "--mutation-probability",
      description = "Default: 1.0 for PAES (mutation-only), 0.2 for the other algorithms.")
  Double mutationProbability;

  @Option(
      names = "--mutation",
      defaultValue = "TM_SAFE_GAP_SHIFT",
      description =
          "Mutation operator: ${COMPLETION-CANDIDATES}. Default: ${DEFAULT-VALUE}. None of them"
              + " breaks a transmembrane helix. RANDOM_CHOICE picks TM_SAFE_GAP_SHIFT or"
              + " MERGE_GAP_RUNS each time.")
  MutationType mutation;

  @Option(
      names = "--initialization",
      defaultValue = "CROSSOVER",
      description =
          "Initial population besides the seeds: ${COMPLETION-CANDIDATES}. Default:"
              + " ${DEFAULT-VALUE} (original). CROSSOVER_DISTINCT avoids duplicates; PERTURBED"
              + " mutates every crossover child.")
  InitializationType initialization;

  @Option(
      names = "--algorithm",
      defaultValue = "NSGA_II",
      description =
          "Search algorithm: ${COMPLETION-CANDIDATES}. Default: ${DEFAULT-VALUE}. NSGA_III needs a"
              + " population size multiple of 4; PAES uses the first seed and ignores the"
              + " crossover; RANDOM_SEARCH ignores the seeds and both operators.")
  AlgorithmType algorithm;

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
                  maxEvaluations,
                  populationSize,
                  populationSize,
                  crossoverProbability,
                  effectiveMutationProbability(),
                  threads,
                  mutation,
                  initialization),
              algorithm,
              mappingPolicy,
              StructuralLddt.DEFAULT_CUTOFF,
              progressInterval);
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
