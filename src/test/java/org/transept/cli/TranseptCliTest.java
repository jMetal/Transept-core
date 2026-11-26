package org.transept.cli;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.transept.io.InvalidInputFileException;
import picocli.CommandLine;

@DisplayName("Unit tests for class TranseptCli")
class TranseptCliTest {

  private StringWriter err;
  private CommandLine commandLine;

  @BeforeEach
  void setUp() {
    err = new StringWriter();
    commandLine = new CommandLine(new TranseptCli()).setErr(new PrintWriter(err));
  }

  @Nested
  @DisplayName("When a command fails")
  class ReportError {

    @Test
    @DisplayName("given invalid input, when reporting it, then one line without stack trace")
    void givenInvalidInput_whenReportingIt_thenOneLineWithoutStackTrace() {
      // Arrange
      var exception = new InvalidInputFileException("Invalid alignment in x.fasta");

      // Act
      int status = TranseptCli.reportError(exception, commandLine, null);

      // Assert
      assertEquals(1, status);
      assertEquals("Error: Invalid alignment in x.fasta" + System.lineSeparator(), err.toString());
    }

    @Test
    @DisplayName("given a file that cannot be read, when reporting it, then its message")
    void givenUnreadableFile_whenReportingIt_thenItsMessage() {
      // Act
      TranseptCli.reportError(new IOException("x.fasta: no such file"), commandLine, null);

      // Assert
      assertEquals("Error: x.fasta: no such file" + System.lineSeparator(), err.toString());
    }

    @Test
    @DisplayName("given an error without message, when reporting it, then its type")
    void givenErrorWithoutMessage_whenReportingIt_thenItsType() {
      // Act
      TranseptCli.reportError(new IllegalArgumentException(), commandLine, null);

      // Assert
      assertEquals("Error: IllegalArgumentException" + System.lineSeparator(), err.toString());
    }

    @Test
    @DisplayName("given an unexpected error, when reporting it, then its stack trace")
    void givenUnexpectedError_whenReportingIt_thenItsStackTrace() {
      // Act
      int status = TranseptCli.reportError(new NullPointerException("x"), commandLine, null);

      // Assert
      assertEquals(1, status);
      assertTrue(err.toString().startsWith("Unexpected error: x"));
      assertTrue(err.toString().contains("\tat "));
    }
  }
}
