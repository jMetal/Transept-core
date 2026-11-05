package org.tmm2structalign.seed;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HexFormat;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.tmm2structalign.alignment.Msa;
import org.tmm2structalign.io.FastaAlignmentReader;
import org.tmm2structalign.io.FastaAlignmentWriter;
import org.tmm2structalign.model.ProteinSequence;
import org.tmm2structalign.model.ProteinSet;

/**
 * Internal seeds computed in advance ({@code seeds} command) and reused by the runs of a dataset.
 *
 * <p>A directory holds the seeds of one dataset ({@code <method>_<dataset>.fasta}) and {@code
 * manifest.tsv}, with the fingerprint of the proteins they align (identifiers, residues and
 * topology) and the SHA-256 of every file. The seeds are used only if the manifest matches the
 * proteins of the run and the files, so seeds of other data are never used silently.
 */
public final class InternalSeedStore {

  /** Name of the manifest file. */
  public static final String MANIFEST = "manifest.tsv";

  /** Directory, under the output root of the runs, of the stored seeds of every dataset. */
  public static final String DIRECTORY = "internal_seeds";

  private static final String FORMAT = "1";

  /**
   * A stored seed.
   *
   * @param method name of the seed, e.g. {@code biojava-g20x2}
   * @param file its FASTA file
   * @param alignment the alignment
   */
  public record StoredSeed(String method, Path file, Msa alignment) {}

  /**
   * What reading a store found.
   *
   * @param seeds the stored seeds, empty if they cannot be used
   * @param problem why they cannot be used, or {@code null} if they can
   */
  public record Lookup(List<StoredSeed> seeds, String problem) {

    /** Returns true if the stored seeds can be used. */
    public boolean usable() {
      return problem == null;
    }
  }

  private InternalSeedStore() {}

  /**
   * Returns the directory of the stored seeds of a dataset.
   *
   * @param outputRoot output root of the runs
   * @param dataset dataset name
   * @return {@code <outputRoot>/internal_seeds/<dataset>}
   */
  public static Path directory(Path outputRoot, String dataset) {
    return outputRoot.resolve(DIRECTORY).resolve(dataset);
  }

  /**
   * Writes the seeds of a dataset and their manifest, replacing earlier ones.
   *
   * @param directory directory of the dataset's seeds
   * @param dataset dataset name
   * @param proteins proteins the seeds align
   * @param seeds the seeds
   * @param generatedBy what generated them, for the manifest (e.g. the commit of the tool)
   * @throws IOException if a file cannot be written
   */
  public static void write(
      Path directory, String dataset, ProteinSet proteins, List<InternalSeeds.Seed> seeds,
      String generatedBy) throws IOException {
    Files.createDirectories(directory);
    StringBuilder manifest = new StringBuilder();
    manifest.append("format\t").append(FORMAT).append('\n');
    manifest.append("dataset\t").append(dataset).append('\n');
    manifest.append("proteins\t").append(fingerprint(proteins)).append('\n');
    manifest.append("generatedBy\t").append(generatedBy).append('\n');
    FastaAlignmentWriter writer = new FastaAlignmentWriter();
    for (InternalSeeds.Seed seed : seeds) {
      String name = seed.method() + "_" + dataset + ".fasta";
      writer.write(seed.alignment(), directory.resolve(name));
      manifest.append("seed\t").append(seed.method()).append('\t').append(name).append('\t')
          .append(sha256(directory.resolve(name))).append('\n');
    }
    Files.writeString(directory.resolve(MANIFEST), manifest, StandardCharsets.UTF_8);
  }

  /**
   * Reads the stored seeds of a dataset, if they match its proteins.
   *
   * @param directory directory of the dataset's seeds
   * @param proteins proteins of the run
   * @return the seeds, or the reason why they cannot be used
   * @throws IOException if a file cannot be read
   */
  public static Lookup read(Path directory, ProteinSet proteins) throws IOException {
    Path manifestFile = directory.resolve(MANIFEST);
    if (!Files.isRegularFile(manifestFile)) {
      return new Lookup(List.of(), "no precomputed internal seeds in " + directory);
    }
    Map<String, String> fields = new LinkedHashMap<>();
    List<String[]> entries = new ArrayList<>();
    for (String line : Files.readAllLines(manifestFile, StandardCharsets.UTF_8)) {
      String[] parts = line.split("\t");
      if (parts.length == 4 && parts[0].equals("seed")) {
        entries.add(parts);
      } else if (parts.length == 2) {
        fields.put(parts[0], parts[1]);
      }
    }
    if (!FORMAT.equals(fields.get("format"))) {
      return new Lookup(List.of(), manifestFile + " has an unknown format");
    }
    if (!fingerprint(proteins).equals(fields.get("proteins"))) {
      return new Lookup(
          List.of(), "the precomputed internal seeds of " + directory + " align other proteins");
    }
    if (entries.isEmpty()) {
      return new Lookup(List.of(), manifestFile + " lists no seeds");
    }
    FastaAlignmentReader reader = new FastaAlignmentReader();
    List<StoredSeed> seeds = new ArrayList<>();
    for (String[] entry : entries) {
      Path file = directory.resolve(entry[2]);
      if (!Files.isRegularFile(file) || !sha256(file).equals(entry[3])) {
        return new Lookup(List.of(), file + " is missing or was changed after it was written");
      }
      seeds.add(new StoredSeed(entry[1], file, reader.read(file, proteins)));
    }
    return new Lookup(List.copyOf(seeds), null);
  }

  /**
   * Fingerprint of a set of proteins: the SHA-256 of their identifiers, residues and topology, in
   * order.
   *
   * @param proteins the proteins
   * @return lower-case hexadecimal digest
   */
  public static String fingerprint(ProteinSet proteins) {
    MessageDigest digest = sha256Digest();
    for (ProteinSequence protein : proteins.proteins()) {
      String line = protein.id() + "\t" + protein.residues() + "\t" + protein.topology() + "\n";
      digest.update(line.getBytes(StandardCharsets.UTF_8));
    }
    return HexFormat.of().formatHex(digest.digest());
  }

  private static String sha256(Path file) throws IOException {
    MessageDigest digest = sha256Digest();
    try (InputStream input = Files.newInputStream(file)) {
      byte[] buffer = new byte[8192];
      for (int read = input.read(buffer); read >= 0; read = input.read(buffer)) {
        digest.update(buffer, 0, read);
      }
    }
    return HexFormat.of().formatHex(digest.digest());
  }

  private static MessageDigest sha256Digest() {
    try {
      return MessageDigest.getInstance("SHA-256");
    } catch (NoSuchAlgorithmException exception) {
      throw new IllegalStateException("SHA-256 is not available", exception);
    }
  }
}
