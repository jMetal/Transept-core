package org.tmm2structalign.objective;

import java.util.List;
import org.tmm2structalign.model.ProteinSet;
import org.tmm2structalign.structure.ProteinStructures;

/** Formulations of the alignment problem: which objectives are optimised, in output order. */
public enum ObjectiveSet {
  /** The lDDT alone, for single-objective algorithms. */
  LDDT,
  /** SOPwTP and lDDT, the published formulation (default). */
  SOPWTP_LDDT,
  /** SOPwTP, lDDT and the aligned-pair coverage. */
  SOPWTP_LDDT_COVERAGE;

  /**
   * Creates the objectives of the formulation.
   *
   * @param proteins proteins of the dataset
   * @param structures structures of the proteins, used by the lDDT
   * @param lddtCutoff inclusion radius of the lDDT, in Å
   * @return the objectives, in output order
   */
  public List<Objective> create(
      ProteinSet proteins, ProteinStructures structures, float lddtCutoff) {
    StructuralLddt lddt = new StructuralLddt(structures, lddtCutoff);
    return switch (this) {
      case LDDT -> List.of(lddt);
      case SOPWTP_LDDT -> List.of(TopologyAwareSumOfPairs.withDefaultParameters(), lddt);
      case SOPWTP_LDDT_COVERAGE ->
          List.of(
              TopologyAwareSumOfPairs.withDefaultParameters(),
              lddt,
              new AlignedPairCoverage(proteins));
    };
  }
}
