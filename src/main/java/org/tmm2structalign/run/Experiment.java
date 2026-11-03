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
import org.tmm2structalign.alignment.Msa;
import org.tmm2structalign.io.FastaAlignmentReader;
import org.tmm2structalign.io.InvalidInputFileException;
import org.tmm2structalign.io.ThreeLineTopologyReader;
import org.tmm2structalign.model.ProteinSet;
import org.tmm2structalign.objective.StructuralLddt;
import org.tmm2structalign.objective.TopologyAwareSumOfPairs;
import org.tmm2structalign.problem.StructuralMsaProblem;
import org.tmm2structalign.structure.ProteinStructures;
import org.tmm2structalign.structure.StructureLoader;

/**
 * Input data of a run, loaded from a {@link RunConfiguration}: proteins, structures, problem and
 * seed alignments.
 *
 * @param problem problem with the SOPwTP and lDDT objectives
 * @param structures structures and mapping report of every protein
 * @param seeds seed alignments, in the configured order
 */
public record Experiment(
    StructuralMsaProblem problem, ProteinStructures structures, List<SeedAlignment> seeds) {

  /**
   * A seed alignment and its provenance.
   *
   * @param method tool that produced the alignment
   * @param file path of the FASTA file
   * @param sha256 hexadecimal SHA-256 of the file content
   * @param alignment the alignment
   */
  public record SeedAlignment(String method, Path file, String sha256, Msa alignment) {}

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
            List.of(
                TopologyAwareSumOfPairs.withDefaultParameters(),
                new StructuralLddt(
                    structures, configuration.lddtCutoff(), configuration.lddtAggregation())));
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
      Path file =
          PrecomputedAlignments.find(
                  configuration.precomputedDirectory(), configuration.dataset(), method)
              .orElseThrow(
                  () ->
                      new InvalidInputFileException(
                          "No unique seed alignment of '" + method + "' in "
                              + configuration.precomputedDirectory()));
      seeds.add(new SeedAlignment(method, file, sha256(file), reader.read(file, proteins)));
    }
    return List.copyOf(seeds);
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
