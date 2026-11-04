package org.tmm2structalign.run;

import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.tmm2structalign.algorithm.AlgorithmConfiguration;
import org.tmm2structalign.algorithm.AlgorithmType;
import org.tmm2structalign.objective.ObjectiveSet;
import org.tmm2structalign.structure.MappingPolicy;

/**
 * Everything needed to reproduce one run.
 *
 * <p>Input paths follow the layout of the {@code Dataset_Structural_GPCRs} repository under
 * {@code dataRoot} (the {@code GPCRdb} directory).
 *
 * @param dataRoot root of the dataset ({@code .../GPCRdb})
 * @param dataset dataset name, e.g. {@code classF}
 * @param seedMethods tools whose precomputed alignments seed the population, in order
 * @param runDirectory output directory of the run; must not exist
 * @param randomSeed seed of the global jMetal generator
 * @param algorithm algorithm parameters
 * @param algorithmType search algorithm
 * @param mappingPolicy policy for mapping sequences onto structures
 * @param lddtCutoff inclusion radius of the structural objective, in Å
 * @param progressInterval evaluations between progress records
 * @param objectives objectives that are optimised
 */
public record RunConfiguration(
    Path dataRoot,
    String dataset,
    List<String> seedMethods,
    Path runDirectory,
    long randomSeed,
    AlgorithmConfiguration algorithm,
    AlgorithmType algorithmType,
    MappingPolicy mappingPolicy,
    float lddtCutoff,
    int progressInterval,
    ObjectiveSet objectives) {

  /** Validates the configuration and freezes the seed list. */
  public RunConfiguration {
    if (dataRoot == null || runDirectory == null || algorithm == null) {
      throw new IllegalArgumentException("dataRoot, runDirectory and algorithm are required");
    }
    if (dataset == null || dataset.isBlank()) {
      throw new IllegalArgumentException("Dataset name cannot be blank");
    }
    if (seedMethods == null || seedMethods.isEmpty()) {
      throw new IllegalArgumentException("At least one seed method is required");
    }
    if (algorithmType == null || mappingPolicy == null || objectives == null) {
      throw new IllegalArgumentException(
          "algorithmType, mappingPolicy and objectives are required");
    }
    if (!Float.isFinite(lddtCutoff) || lddtCutoff <= 0f) {
      throw new IllegalArgumentException("lddtCutoff must be positive: " + lddtCutoff);
    }
    if (progressInterval <= 0) {
      throw new IllegalArgumentException("progressInterval must be positive: " + progressInterval);
    }
    seedMethods = List.copyOf(seedMethods);
  }

  /**
   * Creates a configuration with the published objectives, SOPwTP and lDDT.
   *
   * @param dataRoot root of the dataset
   * @param dataset dataset name
   * @param seedMethods tools whose precomputed alignments seed the population
   * @param runDirectory output directory of the run; must not exist
   * @param randomSeed seed of the global jMetal generator
   * @param algorithm algorithm parameters
   * @param algorithmType search algorithm
   * @param mappingPolicy policy for mapping sequences onto structures
   * @param lddtCutoff inclusion radius of the structural objective, in Å
   * @param progressInterval evaluations between progress records
   */
  public RunConfiguration(
      Path dataRoot,
      String dataset,
      List<String> seedMethods,
      Path runDirectory,
      long randomSeed,
      AlgorithmConfiguration algorithm,
      AlgorithmType algorithmType,
      MappingPolicy mappingPolicy,
      float lddtCutoff,
      int progressInterval) {
    this(
        dataRoot,
        dataset,
        seedMethods,
        runDirectory,
        randomSeed,
        algorithm,
        algorithmType,
        mappingPolicy,
        lddtCutoff,
        progressInterval,
        ObjectiveSet.SOPWTP_LDDT);
  }

  /**
   * Creates a configuration that uses NSGA-II, the published algorithm.
   *
   * @param dataRoot root of the dataset
   * @param dataset dataset name
   * @param seedMethods tools whose precomputed alignments seed the population
   * @param runDirectory output directory of the run; must not exist
   * @param randomSeed seed of the global jMetal generator
   * @param algorithm algorithm parameters
   * @param mappingPolicy policy for mapping sequences onto structures
   * @param lddtCutoff inclusion radius of the structural objective, in Å
   * @param progressInterval evaluations between progress records
   */
  public RunConfiguration(
      Path dataRoot,
      String dataset,
      List<String> seedMethods,
      Path runDirectory,
      long randomSeed,
      AlgorithmConfiguration algorithm,
      MappingPolicy mappingPolicy,
      float lddtCutoff,
      int progressInterval) {
    this(
        dataRoot,
        dataset,
        seedMethods,
        runDirectory,
        randomSeed,
        algorithm,
        AlgorithmType.NSGA_II,
        mappingPolicy,
        lddtCutoff,
        progressInterval);
  }

  /** Returns the 3-line topology file of the dataset. */
  public Path topologyFile() {
    return dataRoot.resolve("sequences").resolve("tmregions")
        .resolve(dataset + "_predicted_topologies.3line");
  }

  /** Returns the directory with the distance matrices of the dataset. */
  public Path distanceDirectory() {
    return dataRoot.resolve("distances").resolve(dataset);
  }

  /** Returns the directory with the AlphaFold PDB files of the dataset. */
  public Path pdbDirectory() {
    return dataRoot.resolve("alphafold_pdb").resolve(dataset);
  }

  /** Returns the directory with the precomputed alignments of the dataset. */
  public Path precomputedDirectory() {
    return dataRoot.resolve("precomputed").resolve(dataset);
  }

  /** Returns the configuration as an ordered map for JSON output. */
  public Map<String, Object> toMap() {
    Map<String, Object> map = new LinkedHashMap<>();
    map.put("dataRoot", dataRoot.toAbsolutePath().toString());
    map.put("dataset", dataset);
    map.put("seedMethods", seedMethods);
    map.put("runDirectory", runDirectory.toAbsolutePath().toString());
    map.put("randomSeed", randomSeed);
    map.put("algorithmType", algorithmType);
    map.put("objectives", objectives);
    map.put("maxEvaluations", algorithm.maxEvaluations());
    map.put("populationSize", algorithm.populationSize());
    map.put("offspringPopulationSize", algorithm.offspringPopulationSize());
    map.put("crossoverProbability", algorithm.crossoverProbability());
    map.put("mutationProbability", algorithm.mutationProbability());
    map.put("mutationType", algorithm.mutationType());
    map.put("threads", algorithm.threads());
    map.put("mappingPolicy", mappingPolicy);
    map.put("lddtCutoff", (double) lddtCutoff);
    map.put("progressInterval", progressInterval);
    return map;
  }
}
