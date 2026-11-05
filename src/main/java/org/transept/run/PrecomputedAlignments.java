package org.transept.run;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Locale;
import java.util.Optional;
import java.util.stream.Stream;

/**
 * Naming rules of the precomputed alignments in {@code precomputed/<dataset>/}.
 *
 * <p>Files are named {@code <tool>_<dataset>.fasta}, but some carry a suffix after the dataset
 * name (e.g. {@code usalign_classF_clean.fasta}). The method name removes the dataset token
 * wherever it appears, so that file is method {@code usalign_clean}.
 */
public final class PrecomputedAlignments {

  private static final String EXTENSION = ".fasta";

  private PrecomputedAlignments() {}

  /**
   * Derives the method name of a precomputed alignment file.
   *
   * @param file FASTA file
   * @param dataset dataset name
   * @return the method name
   */
  public static String methodName(Path file, String dataset) {
    String stem = file.getFileName().toString();
    if (stem.toLowerCase(Locale.ROOT).endsWith(EXTENSION)) {
      stem = stem.substring(0, stem.length() - EXTENSION.length());
    }
    String token = "_" + dataset;
    int position = stem.toLowerCase(Locale.ROOT).indexOf(token.toLowerCase(Locale.ROOT));
    return position < 0
        ? stem
        : stem.substring(0, position) + stem.substring(position + token.length());
  }

  /**
   * Lists the FASTA files of a directory, sorted by name.
   *
   * @param directory directory to list
   * @return FASTA files
   * @throws IOException if the directory cannot be listed
   */
  public static List<Path> list(Path directory) throws IOException {
    try (Stream<Path> listing = Files.list(directory)) {
      return listing
          .filter(
              path -> path.getFileName().toString().toLowerCase(Locale.ROOT).endsWith(EXTENSION))
          .sorted()
          .toList();
    }
  }

  /**
   * Finds the file of a method.
   *
   * @param directory {@code precomputed/<dataset>/}
   * @param dataset dataset name
   * @param method method name as returned by {@link #methodName(Path, String)}
   * @return the file, if exactly one matches
   * @throws IOException if the directory cannot be listed
   */
  public static Optional<Path> find(Path directory, String dataset, String method)
      throws IOException {
    Optional<Path> found = Optional.empty();
    if (Files.isDirectory(directory)) {
      List<Path> matches =
          list(directory).stream()
              .filter(file -> methodName(file, dataset).equals(method))
              .toList();
      found = matches.size() == 1 ? Optional.of(matches.get(0)) : Optional.empty();
    }
    return found;
  }
}
