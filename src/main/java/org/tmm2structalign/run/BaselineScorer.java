package org.tmm2structalign.run;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import org.tmm2structalign.alignment.GapColumnCompactor;
import org.tmm2structalign.alignment.Msa;
import org.tmm2structalign.io.FastaAlignmentReader;
import org.tmm2structalign.io.InvalidInputFileException;
import org.tmm2structalign.io.ThreeLineTopologyReader;
import org.tmm2structalign.model.ProteinSet;
import org.tmm2structalign.objective.Objective;
import org.tmm2structalign.objective.StructuralLddt;
import org.tmm2structalign.objective.TopologyAwareSumOfPairs;
import org.tmm2structalign.structure.MappingPolicy;
import org.tmm2structalign.structure.ProteinStructures;
import org.tmm2structalign.structure.StructureLoader;

/**
 * Scores precomputed alignments with the project objectives, to compare tools with the fronts
 * found by the algorithm. Replaces {@code BaselineMain}: the dataset root is a parameter and
 * output paths are relative to it.
 */
public final class BaselineScorer {

  /**
   * Score of one alignment file.
   *
   * @param dataset dataset name
   * @param method tool name derived from the file name
   * @param file path of the file relative to the dataset root
   * @param values objective values in their natural direction; empty on failure
   * @param error failure message; empty on success
   */
  public record BaselineScore(
      String dataset, String method, String file, List<Double> values, String error) {

    /** Returns true if the alignment could be scored. */
    public boolean isScored() {
      return error.isEmpty();
    }
  }

  private final Path dataRoot;
  private final MappingPolicy mappingPolicy;
  private final float cutoff;

  /**
   * Creates a scorer.
   *
   * @param dataRoot dataset root ({@code .../GPCRdb})
   * @param mappingPolicy sequence-to-structure mapping policy
   * @param cutoff inclusion radius of the structural objective, in Å
   */
  public BaselineScorer(
      Path dataRoot, MappingPolicy mappingPolicy, float cutoff) {
    this.dataRoot = dataRoot;
    this.mappingPolicy = mappingPolicy;
    this.cutoff = cutoff;
  }

  /**
   * Returns the names of the objectives, in the order of {@link BaselineScore#values()}.
   *
   * @return objective names
   */
  public static List<String> objectiveNames() {
    return List.of(TopologyAwareSumOfPairs.NAME, StructuralLddt.NAME);
  }

  /**
   * Scores every {@code *.fasta} file of {@code precomputed/<dataset>/}, sorted by name.
   *
   * @param dataset dataset name
   * @return one score per file; unreadable or inconsistent files are reported, not thrown
   * @throws IOException if the topology, structures or directory cannot be read
   */
  public List<BaselineScore> score(String dataset) throws IOException {
    Path directory = dataRoot.resolve("precomputed").resolve(dataset);
    if (!Files.isDirectory(directory)) {
      throw new InvalidInputFileException("Missing precomputed directory: " + directory);
    }
    Path topology =
        dataRoot.resolve("sequences/tmregions").resolve(dataset + "_predicted_topologies.3line");
    ProteinSet proteins = new ThreeLineTopologyReader().read(topology);
    Path pdbDirectory = dataRoot.resolve("alphafold_pdb").resolve(dataset);
    ProteinStructures structures =
        new StructureLoader(mappingPolicy)
            .load(
                proteins,
                dataRoot.resolve("distances").resolve(dataset),
                Files.isDirectory(pdbDirectory) ? Optional.of(pdbDirectory) : Optional.empty());
    List<Objective> objectives =
        List.of(
            TopologyAwareSumOfPairs.withDefaultParameters(),
            new StructuralLddt(structures, cutoff));

    List<Path> files = PrecomputedAlignments.list(directory);
    List<BaselineScore> scores = new ArrayList<>(files.size());
    for (Path file : files) {
      scores.add(scoreFile(dataset, file, proteins, objectives));
    }
    return scores;
  }

  private BaselineScore scoreFile(
      String dataset, Path file, ProteinSet proteins, List<Objective> objectives)
      throws IOException {
    String method = PrecomputedAlignments.methodName(file, dataset);
    String relative = dataRoot.relativize(file).toString();
    BaselineScore score;
    try {
      Msa read = new FastaAlignmentReader().read(file, proteins);
      Msa msa = Msa.of(proteins, GapColumnCompactor.compact(read.rows()));
      List<Double> values = objectives.stream().map(objective -> objective.evaluate(msa)).toList();
      score = new BaselineScore(dataset, method, relative, values, "");
    } catch (InvalidInputFileException exception) {
      score = new BaselineScore(dataset, method, relative, List.of(), exception.getMessage());
    }
    return score;
  }
}
