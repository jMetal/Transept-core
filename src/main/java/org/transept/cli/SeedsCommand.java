package org.transept.cli;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import org.transept.io.ThreeLineTopologyReader;
import org.transept.model.ProteinSet;
import org.transept.run.BuildInfo;
import org.transept.seed.InternalSeedStore;
import org.transept.seed.InternalSeeds;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

/**
 * Computes the internal seeds of datasets in advance, so that the runs that use {@code --seeds
 * internal} with the same output root read them instead of computing them every time.
 */
@Command(
    name = "seeds",
    mixinStandardHelpOptions = true,
    description = "Compute the internal seed alignments of datasets in advance.")
public final class SeedsCommand implements Callable<Integer> {

  @Option(names = "--data-root", required = true, description = "Dataset root (.../GPCRdb).")
  Path dataRoot;

  @Option(
      names = "--datasets",
      required = true,
      split = ",",
      description = "Dataset names, e.g. classF,classT2_19.")
  List<String> datasets;

  @Option(
      names = "--output-root",
      required = true,
      description =
          "Root of the results tree, as in run: the seeds go to <output-root>/internal_seeds/"
              + "<dataset>/.")
  Path outputRoot;

  @Option(
      names = "--count",
      defaultValue = "10",
      description = "Distinct seeds per dataset, at most 10. Default: ${DEFAULT-VALUE}.")
  int count;

  @Override
  public Integer call() throws IOException {
    Map<String, Object> build = BuildInfo.describe();
    String generatedBy = "Transept " + build.get("gitCommit");
    for (String dataset : datasets) {
      Path topology =
          dataRoot.resolve("sequences").resolve("tmregions")
              .resolve(dataset + "_predicted_topologies.3line");
      if (!Files.isRegularFile(topology)) {
        System.err.println("Missing topology file: " + topology);
        return 1;
      }
      ProteinSet proteins = new ThreeLineTopologyReader().read(topology);
      long start = System.currentTimeMillis();
      List<InternalSeeds.Seed> seeds = InternalSeeds.generate(proteins, count);
      Path directory = InternalSeedStore.directory(outputRoot, dataset);
      InternalSeedStore.write(directory, dataset, proteins, seeds, generatedBy);
      System.out.println(
          dataset + ": " + seeds.size() + " seeds in " + directory + " ("
              + (System.currentTimeMillis() - start) / 1000.0 + " s)");
    }
    return 0;
  }
}
