package org.transept.seed;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.AlignedSequence;
import org.biojava.nbio.core.alignment.template.Profile;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.transept.alignment.GapColumnCompactor;
import org.transept.alignment.GapIntervals;
import org.transept.alignment.Msa;
import org.transept.model.ProteinSet;

/**
 * Seed alignments computed with the progressive aligner of BioJava (a guide tree and profile
 * alignment, in the style of ClustalW, with BLOSUM62), so that a run needs no external aligner.
 *
 * <p>BioJava is deterministic, so different seeds come from different gap penalties: the
 * settings of {@link #SETTINGS} are tried in order, alignments equal to an earlier one (once
 * gap-only columns are removed) are skipped, and the generation stops at the requested number of
 * distinct seeds. The settings avoid low gap-opening penalties with long extensions, which gave
 * poor alignments with many gaps inside helices on the GPCR benchmark.
 */
public final class InternalSeeds {

  /** Name of the internal seeds in {@code --seeds}. */
  public static final String NAME = "internal";

  /** Largest number of internal seeds. */
  public static final int MAX_SEEDS = 10;

  /**
   * Gap penalties of the aligner.
   *
   * @param gapOpening penalty for opening a gap
   * @param gapExtension penalty for extending it
   */
  public record Setting(int gapOpening, int gapExtension) {

    /** Returns the name of the seed, e.g. {@code biojava-g20x2}. */
    public String method() {
      return "biojava-g" + gapOpening + "x" + gapExtension;
    }
  }

  /** Gap penalties tried, in order. */
  public static final List<Setting> SETTINGS =
      List.of(
          new Setting(20, 2),
          new Setting(30, 2),
          new Setting(15, 1),
          new Setting(18, 1),
          new Setting(25, 2),
          new Setting(40, 2),
          new Setting(10, 1),
          new Setting(20, 1),
          new Setting(30, 1),
          new Setting(12, 1),
          new Setting(25, 1),
          new Setting(35, 2));

  /**
   * A generated seed.
   *
   * @param method name of the seed (the aligner and its setting)
   * @param alignment the alignment
   */
  public record Seed(String method, Msa alignment) {}

  private InternalSeeds() {}

  /**
   * Generates distinct seed alignments of a set of proteins.
   *
   * @param proteins proteins to align
   * @param count number of distinct seeds wanted; fewer are returned if the settings give fewer
   *     distinct alignments
   * @return the seeds, in the order of {@link #SETTINGS}
   * @throws IllegalArgumentException if a sequence has residues that BioJava does not accept
   */
  public static List<Seed> generate(ProteinSet proteins, int count) {
    if (count < 1 || count > MAX_SEEDS) {
      throw new IllegalArgumentException("count must be in [1, " + MAX_SEEDS + "]: " + count);
    }
    List<org.biojava.nbio.core.sequence.ProteinSequence> sequences = sequences(proteins);
    List<Seed> seeds = new ArrayList<>();
    Set<List<GapIntervals>> seen = new HashSet<>();
    for (Setting setting : SETTINGS) {
      if (seeds.size() == count) {
        break;
      }
      Msa alignment = align(proteins, sequences, setting);
      if (seen.add(GapColumnCompactor.compact(alignment.rows()))) {
        seeds.add(new Seed(setting.method(), alignment));
      }
    }
    return List.copyOf(seeds);
  }

  private static List<org.biojava.nbio.core.sequence.ProteinSequence> sequences(
      ProteinSet proteins) {
    List<org.biojava.nbio.core.sequence.ProteinSequence> sequences = new ArrayList<>();
    for (int index = 0; index < proteins.size(); index++) {
      try {
        var sequence =
            new org.biojava.nbio.core.sequence.ProteinSequence(proteins.get(index).residues());
        sequence.setOriginalHeader(Integer.toString(index));
        sequences.add(sequence);
      } catch (CompoundNotFoundException exception) {
        throw new IllegalArgumentException(
            "BioJava cannot align " + proteins.get(index).id() + ": " + exception.getMessage(),
            exception);
      }
    }
    return sequences;
  }

  private static Msa align(
      ProteinSet proteins,
      List<org.biojava.nbio.core.sequence.ProteinSequence> sequences,
      Setting setting) {
    Profile<org.biojava.nbio.core.sequence.ProteinSequence, AminoAcidCompound> profile =
        Alignments.getMultipleSequenceAlignment(
            sequences,
            new SimpleGapPenalty(setting.gapOpening(), setting.gapExtension()),
            SubstitutionMatrixHelper.getBlosum62());
    String[] rows = new String[proteins.size()];
    for (AlignedSequence<org.biojava.nbio.core.sequence.ProteinSequence, AminoAcidCompound> row :
        profile.getAlignedSequences()) {
      int index = Integer.parseInt(row.getOriginalSequence().getOriginalHeader());
      rows[index] = row.getSequenceAsString();
    }
    return Msa.fromAlignedRows(proteins, List.of(rows));
  }
}
