package org.tmm2structalign.run;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HexFormat;
import java.util.List;
import java.util.Optional;
import java.util.logging.Logger;
import org.tmm2structalign.alignment.Msa;
import org.tmm2structalign.io.FastaAlignmentReader;
import org.tmm2structalign.io.FastaAlignmentWriter;
import org.tmm2structalign.io.InvalidInputFileException;
import org.tmm2structalign.io.ThreeLineTopologyReader;
import org.tmm2structalign.model.ProteinSet;
import org.tmm2structalign.operator.TmHelixRepair;
import org.tmm2structalign.problem.StructuralMsaProblem;
import org.tmm2structalign.seed.InternalSeedStore;
import org.tmm2structalign.seed.InternalSeeds;
import org.tmm2structalign.structure.ProteinStructures;
import org.tmm2structalign.structure.StructureLoader;

/**
 * Input data of a run, loaded from a {@link RunConfiguration}: proteins, structures, problem and
 * seed alignments.
 *
 * @param problem problem with the configured objectives
 * @param structures structures and mapping report of every protein
 * @param seeds seed alignments, in the configured order
 */
public record Experiment(
    StructuralMsaProblem problem, ProteinStructures structures, List<SeedAlignment> seeds) {

  private static final Logger LOGGER = Logger.getLogger(Experiment.class.getName());

  /** Directory of the run where the internal seeds are written. */
  public static final String INTERNAL_SEEDS_DIRECTORY = "seeds";

  /**
   * A seed alignment and its provenance.
   *
   * @param method tool that produced the alignment ({@code internal:<setting>} for the seeds
   *     computed by the tool)
   * @param file path of the FASTA file (for internal seeds, the copy written to the run directory)
   * @param sha256 hexadecimal SHA-256 of the file content
   * @param alignment the alignment, repaired if the configuration asks for it
   * @param repairedGapRuns gap runs moved out of the helices by the repair
   * @param origin {@code file} (an alignment of the dataset), {@code precomputed} (internal seed
   *     computed in advance) or {@code computed} (internal seed computed in the run)
   */
  public record SeedAlignment(
      String method, Path file, String sha256, Msa alignment, int repairedGapRuns,
      String origin) {}

  /**
   * Loads every input of a run.
   *
   * @param configuration run configuration
   * @return the loaded experiment
   * @throws IOException if a file cannot be read
   * @throws InvalidInputFileException if a file is missing or inconsistent; every configured
   *     seed must exist and match the proteins
   */
  public static Experiment load(RunConfiguration configuration) throws IOException {
    requireFile(configuration.topologyFile(), "topology file");
    ProteinSet proteins = new ThreeLineTopologyReader().read(configuration.topologyFile());
    Path pdbDirectory = configuration.pdbDirectory();
    ProteinStructures structures =
        new StructureLoader(configuration.mappingPolicy())
            .load(
                proteins,
                configuration.distanceDirectory(),
                Files.isDirectory(pdbDirectory) ? Optional.of(pdbDirectory) : Optional.empty());
    StructuralMsaProblem problem =
        new StructuralMsaProblem(
            configuration.dataset(),
            proteins,
            configuration
                .objectives()
                .create(proteins, structures, configuration.lddtCutoff()));
    return new Experiment(problem, structures, loadSeeds(configuration, proteins));
  }

  /** Returns the seed alignments without provenance. */
  public List<Msa> seedAlignments() {
    return seeds.stream().map(SeedAlignment::alignment).toList();
  }

  private static List<SeedAlignment> loadSeeds(RunConfiguration configuration, ProteinSet proteins)
      throws IOException {
    FastaAlignmentReader reader = new FastaAlignmentReader();
    List<SeedAlignment> seeds = new ArrayList<>();
    for (String method : configuration.seedMethods()) {
      if (InternalSeeds.NAME.equals(method)) {
        seeds.addAll(internalSeeds(configuration, proteins));
        continue;
      }
      Path file =
          PrecomputedAlignments.find(
                  configuration.precomputedDirectory(), configuration.dataset(), method)
              .orElseThrow(
                  () ->
                      new InvalidInputFileException(
                          "No unique seed alignment of '" + method + "' in "
                              + configuration.precomputedDirectory()));
      seeds.add(seed(configuration, method, file, reader.read(file, proteins), "file"));
    }
    return List.copyOf(seeds);
  }

  /**
   * Reads the precomputed internal seeds if they match the proteins; otherwise computes them and
   * writes them to {@code seeds/} of the run directory.
   */
  private static List<SeedAlignment> internalSeeds(
      RunConfiguration configuration, ProteinSet proteins) throws IOException {
    if (configuration.internalSeedsDirectory() != null) {
      InternalSeedStore.Lookup lookup =
          InternalSeedStore.read(configuration.internalSeedsDirectory(), proteins);
      if (lookup.usable()) {
        List<SeedAlignment> seeds = new ArrayList<>();
        for (InternalSeedStore.StoredSeed stored : lookup.seeds()) {
          seeds.add(
              seed(configuration, InternalSeeds.NAME + ":" + stored.method(), stored.file(),
                  stored.alignment(), "precomputed"));
        }
        return seeds;
      }
      LOGGER.info(() -> "Computing the internal seeds: " + lookup.problem());
    }
    Path directory = configuration.runDirectory().resolve(INTERNAL_SEEDS_DIRECTORY);
    Files.createDirectories(directory);
    List<SeedAlignment> seeds = new ArrayList<>();
    for (InternalSeeds.Seed generated : InternalSeeds.generate(proteins, InternalSeeds.MAX_SEEDS)) {
      Path file = directory.resolve(generated.method() + ".fasta");
      new FastaAlignmentWriter().write(generated.alignment(), file);
      seeds.add(
          seed(configuration, InternalSeeds.NAME + ":" + generated.method(), file,
              generated.alignment(), "computed"));
    }
    return seeds;
  }

  private static SeedAlignment seed(
      RunConfiguration configuration, String method, Path file, Msa alignment, String origin)
      throws IOException {
    if (!configuration.repairSeeds()) {
      return new SeedAlignment(method, file, sha256(file), alignment, 0, origin);
    }
    TmHelixRepair.Result repaired = TmHelixRepair.repair(alignment);
    return new SeedAlignment(
        method, file, sha256(file), repaired.alignment(), repaired.movedRuns(), origin);
  }

  private static void requireFile(Path file, String description) {
    if (!Files.isRegularFile(file)) {
      throw new InvalidInputFileException("Missing " + description + ": " + file);
    }
  }

  /**
   * Computes the SHA-256 digest of a file.
   *
   * @param file file to hash
   * @return lower-case hexadecimal digest
   * @throws IOException if the file cannot be read
   */
  public static String sha256(Path file) throws IOException {
    try (InputStream input = Files.newInputStream(file)) {
      MessageDigest digest = MessageDigest.getInstance("SHA-256");
      byte[] buffer = new byte[8192];
      int read = input.read(buffer);
      while (read >= 0) {
        digest.update(buffer, 0, read);
        read = input.read(buffer);
      }
      return HexFormat.of().formatHex(digest.digest());
    } catch (NoSuchAlgorithmException exception) {
      throw new IllegalStateException("SHA-256 is not available", exception);
    }
  }
}
