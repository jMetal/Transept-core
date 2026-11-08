package org.transept.operator;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.transept.alignment.GapIntervals;
import org.transept.alignment.MsaSolution;
import org.transept.model.ProteinSequence;
import org.transept.model.ProteinSet;
import org.transept.testing.Alignments;
import org.transept.testing.SeededRandom;

@DisplayName("Unit tests for class TmSafeUniformCrossover")
class TmSafeUniformCrossoverTest {

  private static final ProteinSet PROTEINS = Alignments.proteins("ACDE:OOOO", "ACDE:OOOO");
  private static final MsaSolution FIRST = Alignments.solution(PROTEINS, "A-CD-E", "AC-DE-");
  private static final MsaSolution SECOND = Alignments.solution(PROTEINS, "AC-DE-", "A-CD-E");

  @Test
  @DisplayName("given two parents, when finding consistent positions, then they agree on every row")
  void givenTwoParents_whenFindingConsistentPositions_thenTheyAgreeOnEveryRow() {
    // Act
    List<int[]> positions = TmSafeUniformCrossover.consistentPositions(FIRST, SECOND);

    // Assert: after columns 0, 1, 3, 4 and 6 of both parents, the rows hold the same residues
    assertEquals(
        List.of("0,0", "1,1", "3,3", "4,4", "6,6"),
        positions.stream().map(p -> p[0] + "," + p[1]).toList());
  }

  @Test
  @DisplayName("given drawn blocks, when crossing, then each child takes its blocks unpadded")
  void givenDrawnBlocks_whenCrossing_thenEachChildTakesItsBlocksWithoutPadding() {
    // Arrange: cross, then blocks from the second, first, second and second parent
    RandomSource random = mock(RandomSource.class);
    when(random.nextDouble()).thenReturn(0.0, 0.9, 0.1, 0.9, 0.9);
    TmSafeUniformCrossover crossover = new TmSafeUniformCrossover(1.0, random);

    // Act
    List<MsaSolution> children = crossover.execute(List.of(FIRST, SECOND));

    // Assert
    assertEquals(List.of("A-CDE-", "AC-D-E"), Alignments.rows(children.get(0)));
    assertEquals(List.of("AC-D-E", "A-CDE-"), Alignments.rows(children.get(1)));
  }

  @Test
  @DisplayName("given parents that agree only at the ends, when crossing, then children are copies")
  void givenParentsThatAgreeOnlyAtTheEnds_whenCrossing_thenChildrenAreCopies() {
    // Arrange
    ProteinSet proteins = Alignments.proteins("AC:OO", "AC:OO");
    MsaSolution first = Alignments.solution(proteins, "AC-", "-AC");
    MsaSolution second = Alignments.solution(proteins, "-AC", "AC-");
    TmSafeUniformCrossover crossover = new TmSafeUniformCrossover(1.0, new SeededRandom(1));

    // Act
    List<MsaSolution> children = crossover.execute(List.of(first, second));

    // Assert
    assertEquals(Alignments.rows(first), Alignments.rows(children.get(0)));
    assertEquals(Alignments.rows(second), Alignments.rows(children.get(1)));
    assertNotSame(first, children.get(0));
  }

  @Test
  @DisplayName("given a failed draw, when crossing, then children are copies of the parents")
  void givenFailedDraw_whenCrossing_thenChildrenAreCopies() {
    TmSafeUniformCrossover crossover = new TmSafeUniformCrossover(0.0, new SeededRandom(1));
    List<MsaSolution> children = crossover.execute(List.of(FIRST, SECOND));
    assertEquals(Alignments.rows(FIRST), Alignments.rows(children.get(0)));
    assertEquals(Alignments.rows(SECOND), Alignments.rows(children.get(1)));
  }

  @Test
  @DisplayName("given random parents, when crossing, then residues stay and no helix newly breaks")
  void givenRandomParents_whenCrossing_thenResiduesAreKeptAndNoHelixNewlyBreaks() {
    // Arrange
    ProteinSet proteins =
        Alignments.proteins(
            "ACDEFGHIKLMN:OOMMMMOOMMMO", "ACDEFGHIKL:OMMMMOOMMO", "ACDEFGHIKLMNPQ:OOOMMMMOOOMMMO");
    Random generator = new Random(7);
    TmSafeUniformCrossover crossover = new TmSafeUniformCrossover(1.0, new SeededRandom(3));
    int exchanges = 0;

    for (int trial = 0; trial < 300; trial++) {
      MsaSolution first = Alignments.randomSolution(proteins, 6, generator);
      MsaSolution second = Alignments.randomSolution(proteins, 6, generator);

      // Act
      List<MsaSolution> children = crossover.execute(List.of(first, second));

      // Assert
      for (MsaSolution child : children) {
        exchanges += Alignments.rows(child).equals(Alignments.rows(first))
                || Alignments.rows(child).equals(Alignments.rows(second))
            ? 0
            : 1;
        for (int sequence = 0; sequence < proteins.size(); sequence++) {
          ProteinSequence protein = proteins.get(sequence);
          GapIntervals row = child.variables().get(sequence);
          assertEquals(protein.length(), row.residueCount());
          Set<Integer> fromParents = new HashSet<>(breaks(protein, first, sequence));
          fromParents.addAll(breaks(protein, second, sequence));
          assertTrue(fromParents.containsAll(breaks(protein, child, sequence)));
        }
      }
    }
    assertTrue(exchanges > 0, "some children must mix both parents");
  }

  @Test
  @DisplayName("given invalid arguments, when creating or executing, then exceptions are thrown")
  void givenInvalidArguments_whenCreatingOrExecuting_thenExceptionsAreThrown() {
    SeededRandom random = new SeededRandom(1);
    assertThrows(IllegalArgumentException.class, () -> new TmSafeUniformCrossover(1.5, random));
    TmSafeUniformCrossover crossover = new TmSafeUniformCrossover(1.0, random);
    assertThrows(IllegalArgumentException.class, () -> crossover.execute(List.of(FIRST)));
    MsaSolution other =
        Alignments.solution(Alignments.proteins("ACDE:OOOO", "ACDE:OOOO"), "ACDE", "ACDE");
    assertThrows(IllegalArgumentException.class, () -> crossover.execute(List.of(FIRST, other)));
    assertEquals(2, crossover.numberOfRequiredParents());
    assertEquals(2, crossover.numberOfGeneratedChildren());
  }

  /** Residues {@code r} of a row with gaps before them and after {@code r - 1}, both TM. */
  private static Set<Integer> breaks(ProteinSequence protein, MsaSolution solution, int sequence) {
    boolean[] gaps = solution.variables().get(sequence).toMask();
    Set<Integer> found = new HashSet<>();
    int residues = 0;
    boolean gapSinceResidue = false;
    for (boolean gap : gaps) {
      if (gap) {
        gapSinceResidue = residues > 0;
      } else {
        if (gapSinceResidue
            && protein.isTransmembrane(residues - 1)
            && protein.isTransmembrane(residues)) {
          found.add(residues);
        }
        residues++;
        gapSinceResidue = false;
      }
    }
    return found;
  }
}
