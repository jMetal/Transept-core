package org.transept.run;

import java.nio.file.Path;

/** Signals that a run would overwrite the results of a previous run. */
public class RunDirectoryExistsException extends IllegalStateException {

  private static final long serialVersionUID = 1L;

  /**
   * Creates the exception.
   *
   * @param directory the existing run directory
   */
  public RunDirectoryExistsException(Path directory) {
    super("Run directory already exists and will not be overwritten: " + directory);
  }
}
