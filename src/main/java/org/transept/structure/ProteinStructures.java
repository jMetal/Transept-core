package org.transept.structure;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Optional;

/**
 * Structural data of every protein of a problem, in canonical order, together with the mapping
 * report of each protein.
 */
public final class ProteinStructures {

  private final List<StructureMapping> mappings;
  private final List<float[][]> distances;

  /**
   * Creates the collection.
   *
   * @param mappings one mapping per protein, in canonical order
   * @param distances one distance matrix per protein; ignored for unusable mappings
   */
  public ProteinStructures(List<StructureMapping> mappings, List<float[][]> distances) {
    if (mappings.size() != distances.size()) {
      throw new IllegalArgumentException("Mappings and distance matrices must match");
    }
    this.mappings = List.copyOf(mappings);
    this.distances = Collections.unmodifiableList(new ArrayList<>(distances));
  }

  /** Returns the number of proteins. */
  public int size() {
    return mappings.size();
  }

  /** Returns the mapping report of a protein. */
  public StructureMapping mapping(int protein) {
    return mappings.get(protein);
  }

  /** Returns the mapping reports in canonical order. */
  public List<StructureMapping> mappings() {
    return mappings;
  }

  /** Returns the distance matrix of a protein whose structure is usable. */
  public Optional<float[][]> distances(int protein) {
    return mappings.get(protein).status().isUsable()
        ? Optional.ofNullable(distances.get(protein))
        : Optional.empty();
  }

  /**
   * Returns the structural row of a residue.
   *
   * @param protein protein index
   * @param residue 0-based residue index
   * @return the row in the distance matrix, or {@code -1} if the residue has no structure
   */
  public int structureIndex(int protein, int residue) {
    int[] map = mappings.get(protein).residueToStructure();
    return residue < map.length ? map[residue] : -1;
  }

  /** Returns the identifiers of proteins whose structure is not used. */
  public List<String> proteinsWithoutStructure() {
    return mappings.stream()
        .filter(mapping -> !mapping.status().isUsable())
        .map(StructureMapping::proteinId)
        .toList();
  }
}
