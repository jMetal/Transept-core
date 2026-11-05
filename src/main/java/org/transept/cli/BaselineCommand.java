package org.transept.cli;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.Callable;
import org.transept.objective.StructuralLddt;
import org.transept.run.BaselineScorer;
import org.transept.run.BaselineScorer.BaselineScore;
import org.transept.structure.MappingPolicy;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

/** Scores the precomputed alignments of one or more datasets and writes them as CSV. */
@Command(
    name = "baseline",
    mixinStandardHelpOptions = true,
    description = "Score every precomputed alignment of the given datasets.")
public final class BaselineCommand implements Callable<Integer> {

  @Option(names = "--data-root", required = true, description = "Dataset root (.../GPCRdb).")
  Path dataRoot;

  @Option(
      names = "--datasets",
      required = true,
      split = ",",
      description = "Dataset names, e.g. classF,classT2_19.")
  List<String> datasets;

  @Option(names = "--output", description = "CSV file to write; standard output if absent.")
  Optional<Path> output;

  @Option(
      names = "--mapping-policy",
      defaultValue = "VERIFIED",
      description =
          "Sequence-structure mapping: ${COMPLETION-CANDIDATES}. Default: ${DEFAULT-VALUE}.")
  MappingPolicy mappingPolicy;

  @Override
  public Integer call() throws IOException {
    BaselineScorer scorer =
        new BaselineScorer(
            dataRoot, mappingPolicy, StructuralLddt.DEFAULT_CUTOFF);
    List<BaselineScore> scores = new ArrayList<>();
    for (String dataset : datasets) {
      scores.addAll(scorer.score(dataset));
    }
    String csv = toCsv(scores);
    if (output.isPresent()) {
      Files.writeString(output.get(), csv, StandardCharsets.UTF_8);
    } else {
      System.out.print(csv);
    }
    return scores.stream().anyMatch(BaselineScore::isScored) ? 0 : 1;
  }

  /**
   * Formats scores as CSV with the header
   * {@code dataset,method,SOPwTP,lDDT,file,status,message}.
   *
   * @param scores scores to format
   * @return CSV text
   */
  static String toCsv(List<BaselineScore> scores) {
    StringBuilder csv =
        new StringBuilder("dataset,method,")
            .append(String.join(",", BaselineScorer.objectiveNames()))
            .append(",file,status,message\n");
    for (BaselineScore score : scores) {
      csv.append(quote(score.dataset())).append(',').append(quote(score.method()));
      for (int index = 0; index < BaselineScorer.objectiveNames().size(); index++) {
        csv.append(',').append(score.isScored() ? score.values().get(index).toString() : "");
      }
      csv.append(',')
          .append(quote(score.file()))
          .append(',')
          .append(score.isScored() ? "ok" : "failed")
          .append(',')
          .append(quote(score.error()))
          .append('\n');
    }
    return csv.toString();
  }

  private static String quote(String value) {
    return "\"" + value.replace("\"", "\"\"") + "\"";
  }
}
