package org.transept.operator;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.stream.LongStream;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.MethodSource;
import org.junit.jupiter.params.provider.ValueSource;
import org.transept.alignment.GapIntervals;
import org.transept.alignment.Msa;
import org.transept.alignment.MsaSolution;
import org.transept.model.ProteinSet;
import org.transept.objective.TopologyAwareSumOfPairs;
import org.transept.testing.Alignments;
import org.transept.testing.SeededRandom;

@DisplayName("Unit tests for class TmSafeWindowRealignMutation")
class TmSafeWindowRealignMutationTest {

  private static List<String> realign(
      ProteinSet proteins, int sequence, int first, int columns, String... rows) {
    MsaSolution solution = Alignments.solution(proteins, rows);
    GapIntervals row = TmSafeWindowRealignMutation.realign(solution, sequence, first, columns);
    solution.variables().set(sequence, row);
    return Alignments.rows(solution);
  }

  @Test
  @DisplayName("given a misplaced gap, when realigning, then the row matches the others")
  void givenMisplacedGap_whenRealigning_thenRowMatchesTheOthers() {
    // Arrange
    ProteinSet proteins = Alignments.proteins("WWCC:OOOO", "WWCC:OOOO", "WWCC:OOOO");

    // Act
    List<String> rows = realign(proteins, 0, 0, 5, "W-WCC", "WWCC-", "WWCC-");

    // Assert: the other rows do not change
    assertEquals(List.of("WWCC-", "WWCC-", "WWCC-"), rows);
  }

  @Test
  @DisplayName("given a window, when realigning, then nothing changes outside it")
  void givenWindow_whenRealigning_thenNothingChangesOutsideIt() {
    // Arrange
    ProteinSet proteins = Alignments.proteins("AWWCC:OOOOO", "AWWCC:OOOOO");

    // Act: the gap at the start stays, and the one inside the window moves
    List<String> rows = realign(proteins, 0, 2, 5, "-AW-WCC", "A-WWCC-");

    // Assert
    assertEquals("-AWWCC-", rows.get(0));
  }

  @Test
  @DisplayName("given a better score across a helix, when realigning, then the helix stays whole")
  void givenBetterScoreAcrossHelix_whenRealigning_thenHelixStaysWhole() {
    // Arrange: the other row would pull a gap between the two transmembrane residues
    ProteinSet proteins = Alignments.proteins("WWCC:OMMO", "WWCC:OOOO");
    MsaSolution solution = Alignments.solution(proteins, "WWCC-", "WW-CC");

    // Act
    GapIntervals row = TmSafeWindowRealignMutation.realign(solution, 0, 0, 5);

    // Assert
    assertEquals(0, Alignments.helixBreakingRuns(proteins.get(0), row));
    assertEquals(5, row.alignedLength());
  }

  @Test
  @DisplayName("given a window inside a helix, when realigning, then nothing changes")
  void givenWindowInsideHelix_whenRealigning_thenNothingChanges() {
    // Arrange: an inherited gap, with no valid place for it inside the window
    ProteinSet proteins = Alignments.proteins("AAAA:MMMM", "AAAAA:OOOOO");

    // Act
    List<String> rows = realign(proteins, 0, 1, 3, "AA-AA", "AAAAA");

    // Assert
    assertEquals("AA-AA", rows.get(0));
  }

  @ParameterizedTest(name = "{0}")
  @CsvSource(
      delimiter = '|',
      value = {
        // Each expected row would be another if one transmembrane residue were enough for PHAT
        "only the mutated residue is transmembrane | PKKD:OMMO | NAWER:OOOOO | PKKD-",
        "only the other residue is transmembrane   | NIEP:OOOO | YYKDL:OMMMO | NIEP-",
        "both helices, with mixed pairs when moved | REYI:OMMO | PQGVI:OMMMO | -REYI"
      })
  @DisplayName("given a topology, when realigning, then the matrices are those of SOPwTP")
  void givenTopology_whenRealigning_thenMatricesAreThoseOfSopwtp(
      String name, String mutated, String other, String expected) {
    // Arrange
    ProteinSet proteins = Alignments.proteins(mutated, other);
    String residues = proteins.get(0).residues();

    // Act
    List<String> rows = realign(proteins, 0, 0, 5, residues + "-", proteins.get(1).residues());

    // Assert
    assertEquals(expected, rows.get(0), name);
    assertEquals(bestScore(proteins, 5, rows), score(proteins, rows), name);
  }

  @ParameterizedTest
  @MethodSource("seeds")
  @DisplayName("given rows without gaps, when realigning one, then it gets the best SOPwTP")
  void givenRowsWithoutGaps_whenRealigningOne_thenItGetsTheBestSopwtp(long seed) {
    // Arrange: with the other rows full and the ends outside the membrane, the score is exact
    Random random = new Random(seed);
    int columns = 9;
    int gaps = 1 + random.nextInt(4);
    ProteinSet proteins =
        Alignments.proteins(
            randomProtein(random, columns - gaps),
            randomProtein(random, columns),
            randomProtein(random, columns));
    String start = proteins.get(0).residues() + "-".repeat(gaps);

    // Act
    List<String> rows =
        realign(proteins, 0, 0, columns, start, proteins.get(1).residues(),
            proteins.get(2).residues());

    // Assert
    assertEquals(bestScore(proteins, columns, rows), score(proteins, rows));
    assertEquals(
        0, Alignments.helixBreakingRuns(proteins.get(0), GapIntervals.fromAlignedRow(rows.get(0))));
  }

  @ParameterizedTest
  @CsvSource({"X", "B", "Z"})
  @DisplayName("given an ambiguous residue facing a helix, when realigning, then it is scored")
  void givenAmbiguousResidueFacingHelix_whenRealigning_thenItIsScored(String residue) {
    // Arrange: PHAT has no ambiguity codes, and SOPwTP scores these pairs with BLOSUM62
    ProteinSet mutated = Alignments.proteins("AWWCC:OMMOO", "A" + residue + "WCCA:OOOOOO");
    ProteinSet other = Alignments.proteins("A" + residue + "WCC:OOOOO", "AWWCCA:OMMOOO");

    // Act & Assert
    assertDoesNotThrow(() -> realign(mutated, 0, 0, 6, "AWWCC-", "A" + residue + "WCCA"));
    assertDoesNotThrow(() -> realign(other, 0, 0, 6, "A" + residue + "WCC-", "AWWCCA"));
  }

  private static LongStream seeds() {
    return LongStream.rangeClosed(1, 300);
  }

  /** A random protein with random topology, whose ends are outside the membrane. */
  private static String randomProtein(Random random, int length) {
    String letters = "ACDEFGHIKLMNPQRSTVWY";
    StringBuilder residues = new StringBuilder();
    StringBuilder topology = new StringBuilder();
    for (int index = 0; index < length; index++) {
      boolean end = index == 0 || index == length - 1;
      residues.append(letters.charAt(random.nextInt(letters.length())));
      topology.append(!end && random.nextBoolean() ? 'M' : 'O');
    }
    return residues + ":" + topology;
  }

  private static double score(ProteinSet proteins, List<String> rows) {
    return TopologyAwareSumOfPairs.withDefaultParameters()
        .evaluate(Msa.fromAlignedRows(proteins, rows));
  }

  /** Best SOPwTP over every placement of the first row that breaks no helix. */
  private static double bestScore(ProteinSet proteins, int columns, List<String> rows) {
    String residues = proteins.get(0).residues();
    double best = Double.NEGATIVE_INFINITY;
    for (int mask = 0; mask < 1 << columns; mask++) {
      if (Integer.bitCount(mask) == residues.length()) {
        StringBuilder row = new StringBuilder();
        int next = 0;
        for (int column = 0; column < columns; column++) {
          row.append((mask >> column & 1) == 1 ? residues.charAt(next++) : '-');
        }
        GapIntervals gaps = GapIntervals.fromAlignedRow(row);
        if (Alignments.helixBreakingRuns(proteins.get(0), gaps) == 0) {
          List<String> candidate = new ArrayList<>(rows);
          candidate.set(0, row.toString());
          best = Math.max(best, score(proteins, candidate));
        }
      }
    }
    return best;
  }

  @ParameterizedTest
  @ValueSource(longs = {1, 2, 3, 4, 5, 6})
  @DisplayName("given random alignments, when mutating, then it stays valid")
  void givenRandomAlignments_whenMutating_thenItIsValidAndBreaksNoMoreHelices(long seed) {
    // Arrange
    ProteinSet proteins =
        Alignments.proteins(
            "ACDEFGHIKLMNPQRS:OOMMMMMMOOMMMMMO",
            "ACDEFGHIKLMN:OMMMMMOOMMMO",
            "ACDEFGHIK:OOOOOOOOO");
    MsaSolution solution = Alignments.randomSolution(proteins, 6, new Random(seed));
    int length = solution.alignmentLength();
    int before = 0;
    for (int row = 0; row < proteins.size(); row++) {
      before += Alignments.helixBreakingRuns(proteins.get(row), solution.variables().get(row));
    }

    // Act
    new TmSafeWindowRealignMutation(1.0, 8, new SeededRandom(seed)).execute(solution);

    // Assert
    int after = 0;
    for (int row = 0; row < proteins.size(); row++) {
      GapIntervals gaps = solution.variables().get(row);
      after += Alignments.helixBreakingRuns(proteins.get(row), gaps);
      assertEquals(proteins.get(row).length(), gaps.residueCount());
      assertEquals(length, gaps.alignedLength());
    }
    assertTrue(after <= before, before + " -> " + after);
  }

  @Test
  @DisplayName("given invalid arguments, when creating or executing, then exceptions are thrown")
  void givenInvalidArguments_whenCreatingOrExecuting_thenExceptionsAreThrown() {
    // Arrange
    SeededRandom random = new SeededRandom(1);

    // Act & Assert
    assertThrows(
        IllegalArgumentException.class, () -> new TmSafeWindowRealignMutation(1.5, 10, random));
    assertThrows(
        IllegalArgumentException.class, () -> new TmSafeWindowRealignMutation(0.5, 1, random));
    assertThrows(
        IllegalArgumentException.class,
        () -> new TmSafeWindowRealignMutation(0.5, 10, random).execute(null));
  }
}
