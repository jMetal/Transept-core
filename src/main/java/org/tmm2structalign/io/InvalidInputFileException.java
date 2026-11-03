package org.tmm2structalign.io;

/** Signals that an input file exists but its content is malformed or inconsistent. */
public class InvalidInputFileException extends RuntimeException {

  private static final long serialVersionUID = 1L;

  /**
   * Creates the exception with a message describing the problem and the file involved.
   *
   * @param message description including the file and the offending record
   */
  public InvalidInputFileException(String message) {
    super(message);
  }

  /**
   * Creates the exception wrapping a lower-level cause.
   *
   * @param message description including the file and the offending record
   * @param cause the underlying exception
   */
  public InvalidInputFileException(String message, Throwable cause) {
    super(message, cause);
  }
}
