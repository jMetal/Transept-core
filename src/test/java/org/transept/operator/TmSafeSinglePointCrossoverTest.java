package org.transept.operator;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
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
import org.transept.alignment.MsaSolution;
import org.transept.model.ProteinSet;
import org.transept.testing.Alignments;
import org.transept.testing.SeededRandom;

@DisplayName("Unit tests for class TmSafeSinglePointCrossover")
class TmSafeSinglePointCrossoverTest {

  private TmSafeSinglePointCrossover crossover;

  @BeforeEach
  void setUp() {
    crossover = new TmSafeSinglePointCrossover(1.0, new SeededRandom(1));
  }

  @Nested
  @DisplayName("When finding compatible cuts")
  class CompatibleCuts {

    @Test
    @DisplayName("given a gapped helix, when finding cuts, then residue boundaries are checked")
    void givenGappedHelix_whenFindingCuts_thenResidueBoundariesAreChecked() {
      // Arrange: column 3 cuts sequence 0 after residue 2, inside the M-M helix.
      ProteinSet proteins = Alignments.proteins("AAAAAA:OMMMMO", "AAAAAA:OMMMMO");
      MsaSolution first = Alignments.solution(proteins, "AA---AAAA", "---AAAAAA");
      MsaSolution second = Alignments.solution(proteins, "A-A-AA-AA", "A--AAA-AA");

      // Act & Assert
      assertEquals(List.of(7), crossover.compatibleCuts(first));
      assertEquals(List.of(7), crossover.compatibleCuts(second));
    }

    @Test
    @DisplayName("given gap-only and terminal columns, when finding cuts, then they are skipped")
    void givenGapOnlyAndTerminalColumns_whenFindingCuts_thenTheyAreSkipped() {
      // Arrange
      ProteinSet proteins = Alignments.proteins("ACD:OOO", "ACD:OOO");
      MsaSolution solution = Alignments.solution(proteins, "A-CD", "A-CD");

      // Act & Assert
      assertEquals(List.of(2), crossover.compatibleCuts(solution));
    }
  }

  @Nested
  @DisplayName("When recombining parents")
  class Recombination {

    @Test
    @DisplayName("given different suffix starts, when recombining, then rows are padded")
    void givenDifferentSuffixStarts_whenRecombining_thenRowsArePadded() {
      // Arrange
      ProteinSet proteins = Alignments.proteins("ACDE:OOOO", "ACDE:OOOO");
      MsaSolution first = Alignments.solution(proteins, "AC--DE", "A-CD-E");
      MsaSolution second = Alignments.solution(proteins, "A-CDE-", "ACD--E");

      // Act
      MsaSolution child = crossover.recombine(first, second, 1);

      // Assert
      assertEquals(List.of("AC--DE-", "A-CD--E"), Alignments.rows(child));
    }

    @Test
    @DisplayName("given a cut, when asking for its padding, then the gaps it adds are counted")
    void givenCut_whenAskingForPadding_thenGapsItAddsAreCounted() {
      // Arrange
      ProteinSet proteins = Alignments.proteins("ACDE:OOOO", "ACDE:OOOO");
      MsaSolution first = Alignments.solution(proteins, "AC--DE", "A-CD-E");
      MsaSolution second = Alignments.solution(proteins, "A-CDE-", "ACD--E");

      // Act & Assert: the child of this cut has "--" of padding in its first row
      assertEquals(2, TmSafeSinglePointCrossover.padding(first, second, 1));
    }

    @Test
    @DisplayName("given several candidates, when crossing, then the least padding is used")
    void givenSeveralCandidates_whenCrossing_thenCutWithLeastPaddingIsUsed() {
      // Arrange
      ProteinSet proteins = Alignments.proteins("ACDE:OOOO", "ACDE:OOOO");
      MsaSolution first = Alignments.solution(proteins, "AC--DE", "A-CD-E");
      MsaSolution second = Alignments.solution(proteins, "A-CDE-", "ACD--E");
      List<Integer> cuts = crossover.compatibleCuts(first);
      int worst = 0;
      int best = 0;
      for (int index = 0; index < cuts.size(); index++) {
        long padding = TmSafeSinglePointCrossover.padding(first, second, cuts.get(index));
        if (padding > TmSafeSinglePointCrossover.padding(first, second, cuts.get(worst))) {
          worst = index;
        }
        if (padding < TmSafeSinglePointCrossover.padding(first, second, cuts.get(best))) {
          best = index;
        }
      }
      RandomSource random = mock(RandomSource.class);
      when(random.nextDouble()).thenReturn(0.0);
      when(random.nextInt(0, cuts.size() - 1)).thenReturn(worst, best);
      TmSafeSinglePointCrossover operator = new TmSafeSinglePointCrossover(1.0, 2, random);

      // Act
      MsaSolution child = operator.execute(List.of(first, second)).get(0);

      // Assert
      assertTrue(worst != best);
      assertEquals(
          Alignments.rows(crossover.recombine(first, second, cuts.get(best))),
          Alignments.rows(child));
    }

    @ParameterizedTest
    @ValueSource(longs = {1, 2, 3})
    @DisplayName("given one candidate, when crossing, then it behaves as without the parameter")
    void givenOneCandidate_whenCrossing_thenItBehavesAsWithoutTheParameter(long seed) {
      // Arrange
      ProteinSet proteins = Alignments.proteins("ACDE:OOOO", "ACDE:OOOO");
      MsaSolution first = Alignments.solution(proteins, "AC--DE", "A-CD-E");
      MsaSolution second = Alignments.solution(proteins, "A-CDE-", "ACD--E");

      // Act
      List<MsaSolution> before =
          new TmSafeSinglePointCrossover(1.0, new SeededRandom(seed))
              .execute(List.of(first, second));
      List<MsaSolution> after =
          new TmSafeSinglePointCrossover(
                  1.0, TmSafeSinglePointCrossover.DEFAULT_CUT_CANDIDATES, new SeededRandom(seed))
              .execute(List.of(first, second));

      // Assert
      assertEquals(Alignments.rows(before.get(0)), Alignments.rows(after.get(0)));
      assertEquals(Alignments.rows(before.get(1)), Alignments.rows(after.get(1)));
    }

    @Test
    @DisplayName("given no safe cut, when crossing, then children are copies of the parents")
    void givenNoSafeCut_whenCrossing_thenChildrenAreCopies() {
      // Arrange
      ProteinSet proteins = Alignments.proteins("AAAAAA:MMMMMM", "AAAAAA:MMMMMM");
      MsaSolution first = Alignments.solution(proteins, "AAAAAA", "AAAAAA");
      MsaSolution second = Alignments.solution(proteins, "AAAAAA-", "-AAAAAA");

      // Act
      List<MsaSolution> children = crossover.execute(List.of(first, second));

      // Assert
      assertEquals(first.variables(), children.get(0).variables());
      assertEquals(second.variables(), children.get(1).variables());
      assertNotSame(first, children.get(0));
      assertNotSame(second, children.get(1));
    }

    @Test
    @DisplayName("given a failed probability draw, when crossing, then children are copies")
    void givenFailedProbabilityDraw_whenCrossing_thenChildrenAreCopies() {
      // Arrange
      RandomSource random = mock(RandomSource.class);
      when(random.nextDouble()).thenReturn(0.99);
      TmSafeSinglePointCrossover rarely = new TmSafeSinglePointCrossover(0.5, random);
      ProteinSet proteins = Alignments.proteins("ACDE:OOOO", "ACDE:OOOO");
      MsaSolution first = Alignments.solution(proteins, "AC--DE", "A-CD-E");
      MsaSolution second = Alignments.solution(proteins, "A-CDE-", "ACD--E");

      // Act
      List<MsaSolution> children = rarely.execute(List.of(first, second));

      // Assert
      assertEquals(first.variables(), children.get(0).variables());
      assertEquals(second.variables(), children.get(1).variables());
    }

    @ParameterizedTest
    @ValueSource(longs = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10})
    @DisplayName("given random parents, when crossing, then children are valid and TM-safe")
    void givenRandomParents_whenCrossing_thenChildrenAreValidAndTmSafe(long seed) {
      // Arrange
      Random random = new Random(seed);
      ProteinSet proteins =
          Alignments.proteins(
              "ACDEFGHIKLMN:OOMMMMOOMMMO",
              "ACDEFGHIKL:OMMMMOOMMO",
              "ACDEFGHIKLMNPQ:OOOMMMMOOMMMMO");
      TmSafeSinglePointCrossover operator =
          new TmSafeSinglePointCrossover(1.0, new SeededRandom(seed));

      for (int iteration = 0; iteration < 20; iteration++) {
        MsaSolution first = Alignments.randomSolution(proteins, 6, random);
        MsaSolution second = Alignments.randomSolution(proteins, 4, random);

        // Act
        List<MsaSolution> children = operator.execute(List.of(first, second));

        // Assert: the constructor validated rows; every cut column keeps helices whole.
        assertEquals(2, children.size());
        for (int cut : operator.compatibleCuts(first)) {
          MsaSolution child = operator.recombine(first, second, cut);
          assertTrue(child.alignmentLength() >= cut + 1);
          assertEquals(
              first.toMsa().alignedRow(0).substring(0, cut + 1),
              child.toMsa().alignedRow(0).substring(0, cut + 1));
        }
      }
    }
  }

  @Nested
  @DisplayName("When the operator is misused")
  class Misuse {

    @Test
    @DisplayName("given parents of different protein sets, when crossing, then it is rejected")
    void givenParentsOfDifferentProteinSets_whenCrossing_thenIllegalArgumentExceptionIsThrown() {
      // Arrange
      MsaSolution first =
          Alignments.solution(Alignments.proteins("AC:OO", "AC:OO"), "AC", "AC");
      MsaSolution second =
          Alignments.solution(Alignments.proteins("AC:MM", "AC:MM"), "AC", "AC");

      // Act & Assert
      assertThrows(
          IllegalArgumentException.class, () -> crossover.execute(List.of(first, second)));
    }

    @Test
    @DisplayName("given one parent, when crossing, then an exception is thrown")
    void givenOneParent_whenCrossing_thenIllegalArgumentExceptionIsThrown() {
      // Arrange
      MsaSolution parent = Alignments.solution(Alignments.proteins("AC:OO", "AC:OO"), "AC", "AC");

      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> crossover.execute(List.of(parent)));
    }

    @ParameterizedTest
    @ValueSource(doubles = {-0.1, 1.1, Double.NaN})
    @DisplayName("given an invalid probability, when creating, then an exception is thrown")
    void givenInvalidProbability_whenCreating_thenIllegalArgumentExceptionIsThrown(double value) {
      // Act & Assert
      assertThrows(IllegalArgumentException.class, () -> new TmSafeSinglePointCrossover(value));
    }

    @Test
    @DisplayName("given no candidate cuts, when creating, then an exception is thrown")
    void givenNoCandidateCuts_whenCreating_thenIllegalArgumentExceptionIsThrown() {
      // Act & Assert
      assertThrows(
          IllegalArgumentException.class,
          () -> new TmSafeSinglePointCrossover(0.5, 0, new SeededRandom(1)));
    }

    @Test
    @DisplayName("given the operator, when querying its arity, then two parents and children")
    void givenOperator_whenQueryingArity_thenTwoParentsAndTwoChildren() {
      // Act & Assert
      assertEquals(2, crossover.numberOfRequiredParents());
      assertEquals(2, crossover.numberOfGeneratedChildren());
      assertEquals(1.0, crossover.crossoverProbability());
    }
  }
}
