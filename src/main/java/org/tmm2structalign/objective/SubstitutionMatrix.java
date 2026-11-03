package org.tmm2structalign.objective;

import java.util.Arrays;

/**
 * Symmetric amino-acid substitution matrix with constant-time lookup.
 *
 * <p>The tables were generated from the original {@code Phat} and {@code Blosum62} classes.
 * The original BLOSUM62 had two typos in the ambiguity rows (B-N and Z-Q were 4 instead of 3);
 * they are corrected here. The GPCRdb sequences contain only the 20 standard residues, so the
 * correction does not change any objective value.
 */
public final class SubstitutionMatrix {

  /** PHAT transmembrane substitution matrix (Ng, Henikoff and Henikoff, 2000). */
  public static final SubstitutionMatrix PHAT =
      parse(
          "PHAT",
          """
                   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
              A    5  -6  -2  -5   1  -3  -5   1  -3   0  -1  -7  -1  -1  -3   2   0  -4  -3   1
              R   -6   9  -3  -7  -8  -2  -6  -5  -4  -6  -6  -1  -6  -7  -7  -6  -6  -7  -6  -7
              N   -2  -3  11   2  -2   2   0  -1   4  -3  -3  -2  -2  -1  -4   1  -1  -5   2  -3
              D   -5  -7   2  12  -7   0   6  -2  -1  -5  -5  -5  -5  -5  -5  -4  -5  -7  -4  -5
              C    1  -8  -2  -7   7  -5  -7  -2  -7  -3  -2 -10  -2   0  -8   1  -1  -4  -1  -2
              Q   -3  -2   2   0  -5   9   1  -2   2  -3  -3  -1  -1  -2  -3  -1  -3   1   0  -3
              E   -5  -6   0   6  -7   1  12  -3  -1  -5  -5  -4  -5  -5  -5  -3  -5  -7  -2  -5
              G    1  -5  -1  -2  -2  -2  -3   9  -4  -2  -2  -5  -1  -2  -3   1  -1  -5  -3  -2
              H   -3  -4   4  -1  -7   2  -1  -4  11  -5  -4  -5  -4  -2  -6  -2  -4  -3   3  -5
              I    0  -6  -3  -5  -3  -3  -5  -2  -5   5   2  -7   3   0  -4  -2  -1  -4  -3   3
              L   -1  -6  -3  -5  -2  -3  -5  -2  -4   2   4  -7   2   1  -5  -2  -1  -3  -2   1
              K   -7  -1  -2  -5 -10  -1  -4  -5  -5  -7  -7   5  -6  -7  -4  -5  -6  -8  -4  -8
              M   -1  -6  -2  -5  -2  -1  -5  -1  -4   3   2  -6   6   0  -5  -2   0  -4  -2   1
              F   -1  -7  -1  -5   0  -2  -5  -2  -2   0   1  -7   0   6  -5  -2  -2   0   4  -1
              P   -3  -7  -4  -5  -8  -3  -5  -3  -6  -4  -5  -4  -5  -5  13  -3  -4  -6  -5  -4
              S    2  -6   1  -4   1  -1  -3   1  -2  -2  -2  -5  -2  -2  -3   6   1  -5  -2  -2
              T    0  -6  -1  -5  -1  -3  -5  -1  -4  -1  -1  -6   0  -2  -4   1   3  -7  -3   0
              W   -4  -7  -5  -7  -4   1  -7  -5  -3  -4  -3  -8  -4   0  -6  -5  -7  11   1  -4
              Y   -3  -6   2  -4  -1   0  -2  -3   3  -3  -2  -4  -2   4  -5  -2  -3   1  11  -3
              V    1  -7  -3  -5  -2  -3  -5  -2  -5   3   1  -8   1  -1  -4  -2   0  -4  -3   4
          """);

  /** BLOSUM62 matrix with the ambiguity codes B, Z and X. */
  public static final SubstitutionMatrix BLOSUM62 =
      parse(
          "BLOSUM62",
          """
                  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X
              A   4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0
              R  -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1
              N  -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1
              D  -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1
              C   0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2
              Q  -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1
              E  -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1
              G   0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1
              H  -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1
              I  -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1
              L  -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1
              K  -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1
              M  -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1
              F  -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1
              P  -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2
              S   1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0
              T   0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0
              W  -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2
              Y  -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1
              V   0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1
              B  -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1
              Z  -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1
              X   0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1
          """);

  private static final int NOT_IN_ALPHABET = -1;

  private final String name;
  private final int[] indexOfLetter;
  private final int[][] scores;

  private SubstitutionMatrix(String name, String alphabet, int[][] scores) {
    this.name = name;
    this.scores = scores;
    this.indexOfLetter = new int[128];
    Arrays.fill(indexOfLetter, NOT_IN_ALPHABET);
    for (int index = 0; index < alphabet.length(); index++) {
      indexOfLetter[alphabet.charAt(index)] = index;
    }
  }

  /** Returns the name of the matrix. */
  public String name() {
    return name;
  }

  /**
   * Returns the substitution score of two residues.
   *
   * @param first upper-case residue letter
   * @param second upper-case residue letter
   * @return the score
   * @throws IllegalArgumentException if a residue is not in the matrix alphabet
   */
  public int score(char first, char second) {
    return scores[index(first)][index(second)];
  }

  private int index(char residue) {
    int index = residue < indexOfLetter.length ? indexOfLetter[residue] : NOT_IN_ALPHABET;
    if (index == NOT_IN_ALPHABET) {
      throw new IllegalArgumentException(
          "Residue '" + residue + "' is not in the alphabet of " + name);
    }
    return index;
  }

  /**
   * Parses a square table whose first line lists the alphabet and whose rows start with the
   * row letter, e.g. {@code "  A R\nA 4 -1\nR -1 5"}.
   */
  private static SubstitutionMatrix parse(String name, String table) {
    String[] lines = table.strip().split("\n");
    String alphabet = String.join("", lines[0].trim().split("\\s+"));
    int[][] scores = new int[alphabet.length()][alphabet.length()];
    for (int row = 0; row < alphabet.length(); row++) {
      String[] fields = lines[row + 1].trim().split("\\s+");
      if (fields.length != alphabet.length() + 1 || fields[0].charAt(0) != alphabet.charAt(row)) {
        throw new IllegalStateException("Malformed row " + row + " of matrix " + name);
      }
      for (int column = 0; column < alphabet.length(); column++) {
        scores[row][column] = Integer.parseInt(fields[column + 1]);
      }
    }
    return new SubstitutionMatrix(name, alphabet, scores);
  }
}
