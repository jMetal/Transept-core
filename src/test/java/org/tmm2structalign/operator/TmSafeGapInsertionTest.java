package org.tmm2structalign.operator;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

import java.util.List;
import java.util.Random;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.tmm2structalign.alignment.MsaSolution;
import org.tmm2structalign.model.ProteinSet;
import org.tmm2structalign.testing.Alignments;
import org.tmm2structalign.testing.SeededRandom;

@DisplayName("Unit tests for the TM-safe gap insertion mutations")
class TmSafeGapInsertionTest {

  private static final ProteinSet HELICES =
      Alignments.proteins("ACDEFGHIKLMNPQRS:OOMMMMMMOOMMMMMO", "ACDEFGHIKLMN:OMMMMMOOMMMO");

  private RandomSource random;

  @BeforeEach
  void setUp() {
    random = mock(RandomSource.class);
  }

  private static int breakingRuns(ProteinSet proteins, MsaSolution solution) {
    int breaking = 0;
    for (int sequence = 0; sequence < proteins.size(); sequence++) {
      breaking +=
          Alignments.helixBreakingRuns(proteins.get(sequence), solution.variables().get(sequence));
    }
    return breaking;
  }

  @Nested
  @DisplayName("When protecting transmembrane helices")
  class HelixProtection {

    @ParameterizedTest
    @ValueSource(ints = {1, 3})
    @DisplayName("given random alignments, when inserting gaps, then no helix is broken")
    void givenRandomAlignments_whenInsertingGaps_thenNoHelixIsBroken(int gapLength) {
      // Arrange
      Random generator = new Random(gapLength);
      List<TmSafeGapInsertion> mutations =
          List.of(
              new TmSafeInsertGapMutation(1.0, gapLength, new SeededRandom(1)),
              new TmSafeSplitResidueRunMutation(1.0, gapLength, new SeededRandom(2)));

      for (TmSafeGapInsertion mutation : mutations) {
        for (int iteration = 0; iteration < 300; iteration++) {
          MsaSolution solution = Alignments.randomSolution(HELICES, 10, generator);
          final int before = breakingRuns(HELICES, solution);
          final int length = solution.alignmentLength();

          // Act
          mutation.execute(solution);

          // Assert: existing breaks may remain, but no new one appears
          assertEquals(before, breakingRuns(HELICES, solution), "iteration " + iteration);
          assertEquals(length + gapLength, solution.alignmentLength());
          assertEquals(2, solution.toMsa().size());
        }
      }
    }
  }

  @Nested
  @DisplayName("When inserting anywhere")
  class InsertGap {

    @Test
    @DisplayName("given a helix, when choosing the position, then only its ends are candidates")
    void givenHelix_whenChoosingPosition_thenOnlyItsEndsAreCandidates() {
      // Arrange: valid positions of ACDEF (O M M M O) are 0, 1, 4 and 5; choose the third (4)
      when(random.nextDouble()).thenReturn(0.0);
      when(random.nextInt(0, 3)).thenReturn(2);
      ProteinSet proteins = Alignments.proteins("ACDEF:OMMMO", "ACDEF:OMMMO");
      MsaSolution solution = Alignments.solution(proteins, "ACDEF", "ACDEF");

      // Act
      new TmSafeInsertGapMutation(1.0, 1, random).execute(solution);

      // Assert
      assertEquals(List.of("ACDE-F", "ACDE-F"), Alignments.rows(solution));
    }

    @Test
    @DisplayName("given a position next to a run, when inserting, then the run grows")
    void givenPositionNextToRun_whenInserting_thenRunGrows() {
      // Arrange: valid positions of AC-D are 0..4; position 2 extends the run
      when(random.nextDouble()).thenReturn(0.0);
      when(random.nextInt(0, 4)).thenReturn(2);
      ProteinSet proteins = Alignments.proteins("ACD:OOO", "ACD:OOO");
      MsaSolution solution = Alignments.solution(proteins, "AC-D", "AC-D");

      // Act
      new TmSafeInsertGapMutation(1.0, 2, random).execute(solution);

      // Assert
      assertEquals(List.of("AC---D", "AC---D"), Alignments.rows(solution));
      assertEquals(1, solution.variables().get(0).runCount());
    }
  }

  @Nested
  @DisplayName("When splitting a run of residues")
  class SplitResidueRun {

    @Test
    @DisplayName("given residues before the first run, when splitting, then they are candidates")
    void givenResiduesBeforeFirstRun_whenSplitting_thenTheyAreCandidates() {
      // Arrange: positions between two residues of ACD--EF are 1, 2 and 6; choose 1, which the
      // original operator never chose
      when(random.nextDouble()).thenReturn(0.0);
      when(random.nextInt(0, 2)).thenReturn(0);
      ProteinSet proteins = Alignments.proteins("ACDEF:OOOOO", "ACDEF:OOOOO");
      MsaSolution solution = Alignments.solution(proteins, "ACD--EF", "ACD--EF");

      // Act
      new TmSafeSplitResidueRunMutation(1.0, 1, random).execute(solution);

      // Assert: a new run is created
      assertEquals(List.of("A-CD--EF", "A-CD--EF"), Alignments.rows(solution));
      assertEquals(2, solution.variables().get(0).runCount());
    }

    @Test
    @DisplayName("given a helix, when splitting, then only boundaries outside it are candidates")
    void givenHelix_whenSplitting_thenOnlyBoundariesOutsideItAreCandidates() {
      // Arrange: between-residue positions of ACDEF (O M M M O) outside the helix: 1 and 4
      when(random.nextDouble()).thenReturn(0.0);
      when(random.nextInt(0, 1)).thenReturn(1);
      ProteinSet proteins = Alignments.proteins("ACDEF:OMMMO", "ACDEF:OMMMO");
      MsaSolution solution = Alignments.solution(proteins, "ACDEF", "ACDEF");

      // Act
      new TmSafeSplitResidueRunMutation(1.0, 1, random).execute(solution);

      // Assert
      assertEquals(List.of("ACDE-F", "ACDE-F"), Alignments.rows(solution));
    }

    @Test
    @DisplayName("given no valid position between residues, when splitting, then it falls back")
    void givenNoValidPositionBetweenResidues_whenSplitting_thenItFallsBackToAnyValidPosition() {
      // Arrange: AC is a helix, so the only position between residues is invalid; valid
      // positions are 0 and 2
      when(random.nextDouble()).thenReturn(0.0);
      when(random.nextInt(0, 1)).thenReturn(1);
      ProteinSet proteins = Alignments.proteins("AC:MM", "AC:MM");
      MsaSolution solution = Alignments.solution(proteins, "AC", "AC");

      // Act
      new TmSafeSplitResidueRunMutation(1.0, 1, random).execute(solution);

      // Assert
      assertEquals(List.of("AC-", "AC-"), Alignments.rows(solution));
    }
  }

  @Nested
  @DisplayName("When checking the probability and the arguments")
  class Arguments {

    @Test
    @DisplayName("given probability zero, when mutating, then the solution is unchanged")
    void givenProbabilityZero_whenMutating_thenSolutionIsUnchanged() {
      // Arrange
      when(random.nextDouble()).thenReturn(0.5);
      ProteinSet proteins = Alignments.proteins("ACD:OOO", "ACD:OOO");
      MsaSolution solution = Alignments.solution(proteins, "AC-D", "ACD-");

      // Act
      MsaSolution result = new TmSafeInsertGapMutation(0.0, 1, random).execute(solution);

      // Assert
      assertSame(solution, result);
      assertEquals(List.of("AC-D", "ACD-"), Alignments.rows(solution));
    }

    @Test
    @DisplayName("given invalid arguments, when creating, then IllegalArgumentException is thrown")
    void givenInvalidArguments_whenCreating_thenIllegalArgumentExceptionIsThrown() {
      assertThrows(
          IllegalArgumentException.class, () -> new TmSafeInsertGapMutation(1.5, 1, random));
      assertThrows(
          IllegalArgumentException.class, () -> new TmSafeSplitResidueRunMutation(0.5, 0, random));
      TmSafeInsertGapMutation mutation = new TmSafeInsertGapMutation(0.5);
      assertThrows(IllegalArgumentException.class, () -> mutation.execute(null));
      assertEquals(1, mutation.gapLength());
      assertEquals(0.5, mutation.mutationProbability());
    }
  }

  @Test
  @DisplayName("given a solution, when mutating, then the residues of every row are preserved")
  void givenSolution_whenMutating_thenResiduesOfEveryRowArePreserved() {
    // Arrange
    MsaSolution solution = Alignments.randomSolution(HELICES, 5, new Random(7));

    // Act
    new TmSafeSplitResidueRunMutation(1.0, 2, new SeededRandom(3)).execute(solution);

    // Assert
    for (int sequence = 0; sequence < HELICES.size(); sequence++) {
      String row = Alignments.rows(solution).get(sequence);
      assertEquals(HELICES.get(sequence).residues(), row.replace("-", ""));
      assertEquals(solution.alignmentLength(), row.length());
    }
  }
}
