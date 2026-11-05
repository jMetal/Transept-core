package org.transept.io;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import org.transept.alignment.Msa;

/** Writes alignments in FASTA format, wrapping sequences at a fixed line width. */
public final class FastaAlignmentWriter {

  /** Line width used by the original implementation. */
  public static final int DEFAULT_LINE_WIDTH = 60;

  private final int lineWidth;

  /** Creates a writer with {@link #DEFAULT_LINE_WIDTH}. */
  public FastaAlignmentWriter() {
    this(DEFAULT_LINE_WIDTH);
  }

  /**
   * Creates a writer with a custom line width.
   *
   * @param lineWidth maximum residues per line, positive
   */
  public FastaAlignmentWriter(int lineWidth) {
    if (lineWidth <= 0) {
      throw new IllegalArgumentException("Line width must be positive: " + lineWidth);
    }
    this.lineWidth = lineWidth;
  }

  /**
   * Formats an alignment as FASTA text, using the protein identifiers as headers.
   *
   * @param msa alignment to format
   * @return FASTA text ending with a line separator
   */
  public String format(Msa msa) {
    StringBuilder text = new StringBuilder();
    for (int sequence = 0; sequence < msa.size(); sequence++) {
      text.append('>').append(msa.proteins().get(sequence).id()).append('\n');
      String aligned = msa.alignedRow(sequence);
      for (int start = 0; start < aligned.length(); start += lineWidth) {
        text.append(aligned, start, Math.min(start + lineWidth, aligned.length())).append('\n');
      }
    }
    return text.toString();
  }

  /**
   * Writes an alignment to a file, replacing any previous content.
   *
   * @param msa alignment to write
   * @param file destination path
   * @throws IOException if the file cannot be written
   */
  public void write(Msa msa, Path file) throws IOException {
    Files.writeString(file, format(msa), StandardCharsets.UTF_8);
  }
}
