package org.transept.cli;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.concurrent.Callable;
import org.transept.alignment.Msa;
import org.transept.io.FastaAlignmentReader;
import org.transept.io.FastaAlignmentWriter;
import org.transept.io.ThreeLineTopologyReader;
import org.transept.model.ProteinSet;
import org.transept.operator.TmHelixRepair;
import picocli.CommandLine.Command;
import picocli.CommandLine.Model.CommandSpec;
import picocli.CommandLine.Option;
import picocli.CommandLine.ParameterException;
import picocli.CommandLine.Parameters;
import picocli.CommandLine.Spec;

/**
 * Moves the gaps that split a transmembrane helix out of alignments, without any search: the same
 * repair that {@code run --repair-seeds} applies to the seed alignments.
 */
@Command(
    name = "repair",
    mixinStandardHelpOptions = true,
    description = {
      "Move the gaps that split a transmembrane helix out of alignments, as run --repair-seeds "
          + "does with the seeds, and write the repaired alignments.",
      "Each gap run inside a helix moves to the nearest position outside the helices; the "
          + "residues keep their order and the length of the alignment does not change."
    })
public final class RepairCommand implements Callable<Integer> {

  @Spec CommandSpec spec;

  @Option(
      names = "--topology",
      required = true,
      description = "Three-line topology file of the proteins of the alignments.")
  Path topology;

  @Option(
      names = "--output-dir",
      required = true,
      description = "Directory for the repaired alignments, written with the same file names.")
  Path outputDirectory;

  @Parameters(
      arity = "1..*",
      paramLabel = "ALIGNMENT",
      description = "Alignments to repair (FASTA), with the proteins of the topology.")
  List<Path> alignments;

  @Override
  public Integer call() throws IOException {
    if (!Files.isRegularFile(topology)) {
      throw new ParameterException(spec.commandLine(), "--topology is not a file: " + topology);
    }
    for (Path alignment : alignments) {
      if (!Files.isRegularFile(alignment)) {
        throw new ParameterException(spec.commandLine(), "Not a file: " + alignment);
      }
      Path target = outputDirectory.resolve(alignment.getFileName());
      if (Files.exists(target) && Files.isSameFile(target, alignment)) {
        throw new ParameterException(
            spec.commandLine(), "--output-dir would overwrite the input " + alignment);
      }
    }
    ProteinSet proteins = new ThreeLineTopologyReader().read(topology);
    Files.createDirectories(outputDirectory);
    FastaAlignmentReader reader = new FastaAlignmentReader();
    FastaAlignmentWriter writer = new FastaAlignmentWriter();
    for (Path alignment : alignments) {
      Path target = outputDirectory.resolve(alignment.getFileName());
      Msa msa = reader.read(alignment, proteins);
      TmHelixRepair.Result repaired = TmHelixRepair.repair(msa);
      writer.write(repaired.alignment(), target);
      System.out.println(
          alignment.getFileName() + ": " + repaired.movedRuns() + " gap runs moved -> " + target);
    }
    return 0;
  }
}
