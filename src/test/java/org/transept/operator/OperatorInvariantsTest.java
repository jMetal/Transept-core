package org.transept.operator;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.stream.Stream;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;
import org.transept.alignment.GapIntervals;
import org.transept.alignment.MsaSolution;
import org.transept.evolver.MsaOperators;
import org.transept.model.ProteinSequence;
import org.transept.model.ProteinSet;
import org.transept.testing.Alignments;
import org.transept.testing.SeededRandom;
import org.uma.jmetal.operator.crossover.CrossoverOperator;
import org.uma.jmetal.operator.mutation.MutationOperator;

/** What every operator must keep, whatever it does: checked on many random alignments. */
@DisplayName("Invariants of every operator")
class OperatorInvariantsTest {

  private static final int ALIGNMENTS = 300;
  private static final String LETTERS = "ACDEFGHIKLMNPQRSTVWYXBZ";

  private static Stream<String> mutations() {
    return MsaOperators.MUTATIONS.stream();
  }

  private static Stream<String> crossovers() {
    return MsaOperators.CROSSOVERS.stream();
  }

  @ParameterizedTest
  @MethodSource("mutations")
  @DisplayName("given random alignments, when mutating, then it is valid and breaks no helix")
  void givenRandomAlignments_whenMutating_thenItIsValidAndBreaksNoHelix(String name) {
    for (long seed = 1; seed <= ALIGNMENTS; seed++) {
      // Arrange
      Random random = new Random(seed);
      ProteinSet proteins = randomProteins(random);
      MsaSolution solution = Alignments.randomSolution(proteins, 8, random);
      int length = solution.alignmentLength();
      List<Set<Integer>> before = brokenBoundaries(solution);
      MutationOperator<MsaSolution> mutation =
          MsaOperators.mutation(
              name, 1.0, random.nextInt(4), 1 + random.nextInt(4), 4 + random.nextInt(20),
              Alignments.fullyStructured(proteins), new SeededRandom(seed));

      // Act
      MsaSolution mutated = mutation.execute(solution);

      // Assert
      assertEquals(length, mutated.alignmentLength(), name + ", seed " + seed);
      assertValid(mutated, before, name + ", seed " + seed);
    }
  }

  @ParameterizedTest
  @MethodSource("crossovers")
  @DisplayName("given random parents, when crossing, then children are valid and break no helix")
  void givenRandomParents_whenCrossing_thenChildrenAreValidAndBreakNoHelix(String name) {
    for (long seed = 1; seed <= ALIGNMENTS; seed++) {
      // Arrange: the second parent shares the proteins, and sometimes part of the alignment
      Random random = new Random(seed);
      ProteinSet proteins = randomProteins(random);
      MsaSolution first = Alignments.randomSolution(proteins, 8, random);
      MsaSolution second =
          random.nextBoolean() ? Alignments.randomSolution(proteins, 5, random) : first.copy();
      new MergeAdjacentGapRunsMutation(1.0, new SeededRandom(seed)).execute(second);
      List<Set<Integer>> allowed = brokenBoundaries(first);
      List<Set<Integer>> ofSecond = brokenBoundaries(second);
      for (int row = 0; row < allowed.size(); row++) {
        allowed.get(row).addAll(ofSecond.get(row));
      }
      final List<String> parents = new ArrayList<>(Alignments.rows(first));
      parents.addAll(Alignments.rows(second));
      CrossoverOperator<MsaSolution> crossover =
          MsaOperators.crossover(name, 1.0, 1 + random.nextInt(8), new SeededRandom(seed));

      // Act
      List<MsaSolution> children = crossover.execute(List.of(first, second));

      // Assert: and the parents do not change
      assertEquals(2, children.size());
      for (MsaSolution child : children) {
        assertValid(child, allowed, name + ", seed " + seed);
      }
      List<String> after = new ArrayList<>(Alignments.rows(first));
      after.addAll(Alignments.rows(second));
      assertEquals(parents, after, name + ", seed " + seed);
    }
  }

  @ParameterizedTest
  @MethodSource("mutations")
  @DisplayName("given repaired alignments, when mutating many times, then no helix is ever split")
  void givenRepairedAlignments_whenMutatingManyTimes_thenNoHelixIsEverSplit(String name) {
    for (long seed = 1; seed <= 30; seed++) {
      // Arrange
      Random random = new Random(seed);
      ProteinSet proteins = randomProteins(random);
      MsaSolution solution =
          new MsaSolution(
              proteins,
              TmHelixRepair.repair(Alignments.randomSolution(proteins, 8, random).toMsa())
                  .alignment()
                  .rows(),
              2);
      MutationOperator<MsaSolution> mutation =
          MsaOperators.mutation(
              name, 1.0, 0, 3, 15, Alignments.fullyStructured(proteins), new SeededRandom(seed));

      // Act
      for (int step = 0; step < 50; step++) {
        solution = mutation.execute(solution);
      }

      // Assert
      for (Set<Integer> broken : brokenBoundaries(solution)) {
        assertTrue(broken.isEmpty(), name + ", seed " + seed);
      }
    }
  }

  private static ProteinSet randomProteins(Random random) {
    String[] proteins = new String[3 + random.nextInt(3)];
    for (int index = 0; index < proteins.length; index++) {
      int length = 12 + random.nextInt(14);
      StringBuilder residues = new StringBuilder();
      StringBuilder topology = new StringBuilder();
      boolean membrane = random.nextBoolean();
      for (int residue = 0; residue < length; residue++) {
        membrane = random.nextInt(5) == 0 ? !membrane : membrane;
        // The ambiguity codes only exist in the matrix used outside the membrane
        int letters = membrane ? LETTERS.length() - 3 : LETTERS.length();
        residues.append(LETTERS.charAt(random.nextInt(letters)));
        topology.append(membrane ? 'M' : 'O');
      }
      proteins[index] = residues + ":" + topology;
    }
    return Alignments.proteins(proteins);
  }

  /** For every row, the residues before each gap run that splits a helix. */
  private static List<Set<Integer>> brokenBoundaries(MsaSolution solution) {
    List<Set<Integer>> broken = new ArrayList<>();
    for (int row = 0; row < solution.variables().size(); row++) {
      ProteinSequence protein = solution.proteins().get(row);
      GapIntervals gaps = solution.variables().get(row);
      int[] before = TmBoundaries.residuesBefore(gaps.toMask());
      Set<Integer> boundaries = new HashSet<>();
      for (int run = 0; run < gaps.runCount(); run++) {
        if (TmBoundaries.splitsHelix(protein, before[gaps.start(run)])) {
          boundaries.add(before[gaps.start(run)]);
        }
      }
      broken.add(boundaries);
    }
    return broken;
  }

  private static void assertValid(
      MsaSolution solution, List<Set<Integer>> allowed, String message) {
    List<Set<Integer>> broken = brokenBoundaries(solution);
    for (int row = 0; row < solution.variables().size(); row++) {
      GapIntervals gaps = solution.variables().get(row);
      assertEquals(solution.proteins().get(row).length(), gaps.residueCount(), message);
      assertEquals(solution.alignmentLength(), gaps.alignedLength(), message);
      assertTrue(allowed.get(row).containsAll(broken.get(row)), message);
    }
  }
}
