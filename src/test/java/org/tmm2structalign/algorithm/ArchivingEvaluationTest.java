package org.tmm2structalign.algorithm;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertSame;

import java.util.List;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.tmm2structalign.alignment.MsaSolution;
import org.tmm2structalign.model.ProteinSet;
import org.tmm2structalign.objective.TopologyAwareSumOfPairs;
import org.tmm2structalign.problem.StructuralMsaProblem;
import org.tmm2structalign.testing.Alignments;
import org.uma.jmetal.component.catalogue.common.evaluation.impl.SequentialEvaluation;
import org.uma.jmetal.util.archive.impl.NonDominatedSolutionListArchive;

@DisplayName("Unit tests for class ArchivingEvaluation")
class ArchivingEvaluationTest {

  @Test
  @DisplayName("given solutions, when evaluating, then copies of them reach the archive")
  void givenSolutions_whenEvaluating_thenCopiesReachTheArchive() {
    // Arrange: with a single objective only the best solution is non-dominated
    ProteinSet proteins = Alignments.proteins("ACDEF:OOMMM", "ACDE:OOMM");
    StructuralMsaProblem problem =
        new StructuralMsaProblem(
            "test", proteins, List.of(TopologyAwareSumOfPairs.withDefaultParameters()));
    MsaSolution aligned = Alignments.solution(proteins, "ACDEF", "ACDE-");
    MsaSolution shifted = Alignments.solution(proteins, "ACDEF----", "----ACDE-");
    var archive = new NonDominatedSolutionListArchive<MsaSolution>();
    var evaluation = new ArchivingEvaluation(new SequentialEvaluation<>(problem), archive);

    // Act
    List<MsaSolution> evaluated = evaluation.evaluate(List.of(aligned, shifted));

    // Assert
    assertEquals(List.of(aligned, shifted), evaluated);
    assertEquals(2, evaluation.computedEvaluations());
    assertSame(problem, evaluation.problem());
    assertSame(archive, evaluation.archive());
    assertEquals(1, archive.size());
    MsaSolution kept = archive.solutions().get(0);
    assertNotSame(aligned, kept);
    assertEquals(aligned.variables(), kept.variables());
  }
}
