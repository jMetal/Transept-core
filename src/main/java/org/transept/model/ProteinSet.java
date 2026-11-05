package org.transept.model;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;

/** Immutable, ordered collection of proteins with unique identifiers. */
public final class ProteinSet {

  private final List<ProteinSequence> proteins;
  private final Map<String, Integer> indexById;

  /**
   * Creates a protein set preserving the given order.
   *
   * @param proteins at least two proteins with distinct identifiers
   * @throws IllegalArgumentException if fewer than two proteins are given or ids repeat
   */
  public ProteinSet(List<ProteinSequence> proteins) {
    if (proteins == null || proteins.size() < 2) {
      throw new IllegalArgumentException("A protein set needs at least two proteins");
    }
    this.proteins = List.copyOf(proteins);
    this.indexById = new HashMap<>();
    for (int index = 0; index < this.proteins.size(); index++) {
      String id = this.proteins.get(index).id();
      if (indexById.put(id, index) != null) {
        throw new IllegalArgumentException("Duplicated protein identifier: " + id);
      }
    }
  }

  /** Returns the number of proteins. */
  public int size() {
    return proteins.size();
  }

  /** Returns the protein at a 0-based index. */
  public ProteinSequence get(int index) {
    return proteins.get(index);
  }

  /** Returns the proteins in their canonical order. */
  public List<ProteinSequence> proteins() {
    return proteins;
  }

  /** Returns the index of a protein identifier, if present. */
  public Optional<Integer> indexOf(String id) {
    return Optional.ofNullable(indexById.get(id));
  }
}
