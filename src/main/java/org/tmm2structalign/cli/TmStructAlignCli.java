package org.tmm2structalign.cli;

import picocli.CommandLine;
import picocli.CommandLine.Command;

/**
 * Entry point of the command-line interface.
 *
 * <p>Arguments can be read from a file with {@code @file}, one option per line, so that a run is
 * fully described by a versionable file.
 */
@Command(
    name = "tm-m2structalign",
    mixinStandardHelpOptions = true,
    versionProvider = VersionProvider.class,
    description = "Structure-guided multiobjective alignment of transmembrane proteins.",
    subcommands = {RunCommand.class, BaselineCommand.class, AlgorithmsCommand.class})
public final class TmStructAlignCli {

  /**
   * Runs the command line and exits with its status code.
   *
   * @param args command-line arguments
   */
  public static void main(String[] args) {
    System.exit(execute(args));
  }

  /**
   * Runs the command line without exiting the JVM.
   *
   * @param args command-line arguments
   * @return 0 on success, non-zero on failure
   */
  public static int execute(String... args) {
    return new CommandLine(new TmStructAlignCli()).execute(args);
  }
}
