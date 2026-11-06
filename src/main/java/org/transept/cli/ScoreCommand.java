package org.transept.cli;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Optional;
import java.util.concurrent.Callable;
import org.transept.evaluation.ScoreReport;
import picocli.CommandLine.Command;
import picocli.CommandLine.Model.CommandSpec;
import picocli.CommandLine.Option;
import picocli.CommandLine.ParameterException;
import picocli.CommandLine.Spec;

/** Compares the alignments of a directory tree with a reference alignment. */
@Command(
    name = "score",
    mixinStandardHelpOptions = true,
    description = "Compare every alignment of a directory with a reference alignment.")
public final class ScoreCommand implements Callable<Integer> {

  @Spec CommandSpec spec;

  @Option(names = "--reference", required = true, description = "Reference alignment (FASTA).")
  Path reference;

  @Option(
      names = "--topology",
      required = true,
      description = "Three-line topology file of the dataset.")
  Path topology;

  @Option(
      names = "--msa-root",
      required = true,
      description = "Directory searched recursively for alignments.")
  Path msaRoot;

  @Option(
      names = "--glob",
      defaultValue = "MSA*.fasta",
      description = "Pattern of the names of the alignment files. Default: ${DEFAULT-VALUE}.")
  String glob;

  @Option(
      names = "--output",
      description = "CSV report to write. Default: evaluation_metrics.csv in the MSA root.")
  Optional<Path> output;

  @Override
  public Integer call() throws IOException {
    requireFile(reference, "--reference");
    requireFile(topology, "--topology");
    if (!Files.isDirectory(msaRoot)) {
      throw new ParameterException(
          spec.commandLine(), "--msa-root is not a directory: " + msaRoot);
    }
    if (glob.isBlank()) {
      throw new ParameterException(spec.commandLine(), "--glob must not be blank");
    }
    Path report = output.orElse(msaRoot.resolve("evaluation_metrics.csv"));
    ScoreReport.Summary summary = ScoreReport.write(reference, topology, msaRoot, glob, report);
    System.out.println(
        "Compared " + summary.scored() + " alignments, " + summary.failed() + " failed: " + report);
    return summary.scored() > 0 ? 0 : 1;
  }

  private void requireFile(Path path, String option) {
    if (!Files.isRegularFile(path)) {
      throw new ParameterException(spec.commandLine(), option + " is not a file: " + path);
    }
  }
}
