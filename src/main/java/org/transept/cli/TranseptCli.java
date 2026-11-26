package org.transept.cli;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.UncheckedIOException;
import org.transept.io.InvalidInputFileException;
import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.ParseResult;

/**
 * Entry point of the command-line interface.
 *
 * <p>Arguments can be read from a file with {@code @file}, one option per line, so that a run is
 * fully described by a versionable file.
 */
@Command(
    name = "transept-core",
    mixinStandardHelpOptions = true,
    versionProvider = VersionProvider.class,
    description = "Structure-guided multiobjective alignment of transmembrane proteins.",
    subcommands = {
      RunCommand.class,
      BaselineCommand.class,
      AlgorithmsCommand.class,
      SeedsCommand.class,
      ScoreCommand.class,
      RepairCommand.class
    })
public final class TranseptCli {

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
    return new CommandLine(new TranseptCli())
        .setExecutionExceptionHandler(TranseptCli::reportError)
        .execute(args);
  }

  /**
   * Reports an error of a command: one line for invalid input, which the user can correct, and
   * the stack trace for anything else, which is a bug.
   */
  static int reportError(Exception exception, CommandLine commandLine, ParseResult parseResult) {
    PrintWriter err = commandLine.getErr();
    if (isInvalidInput(exception)) {
      err.println("Error: " + describe(exception));
    } else {
      err.println("Unexpected error: " + describe(exception));
      exception.printStackTrace(err);
    }
    err.flush();
    return commandLine.getCommandSpec().exitCodeOnExecutionException();
  }

  private static boolean isInvalidInput(Exception exception) {
    return exception instanceof IllegalArgumentException
        || exception instanceof IllegalStateException
        || exception instanceof InvalidInputFileException
        || exception instanceof IOException
        || exception instanceof UncheckedIOException;
  }

  private static String describe(Exception exception) {
    String message = exception.getMessage();
    return message == null || message.isBlank() ? exception.getClass().getSimpleName() : message;
  }
}
