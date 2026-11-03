package org.tmm2structalign.structure;

/**
 * Immutable sparse representation of the local structural neighbourhood of a protein: for every
 * structural residue, only the neighbours closer than a cutoff are stored, with their distance
 * and a distance-decay weight.
 *
 * <p>Port of {@code org.tm_msaligner.structure.StructuralNeighborhood} plus the construction
 * code of {@code LDDTStructuralScore}; arithmetic is kept in {@code float} as in the original.
 */
public final class StructuralNeighborhood {

  private static final float MIN_CONTACT_WEIGHT = 1.0e-6f;

  private final int residueCount;
  private final int[][] neighborIndices;
  private final float[][] neighborDistances;
  private final float[][] neighborWeights;

  private StructuralNeighborhood(
      int[][] neighborIndices, float[][] neighborDistances, float[][] neighborWeights) {
    this.residueCount = neighborIndices.length;
    this.neighborIndices = neighborIndices;
    this.neighborDistances = neighborDistances;
    this.neighborWeights = neighborWeights;
  }

  /**
   * Builds the neighbourhood of every residue from a dense distance matrix.
   *
   * @param distances square Cα distance matrix in Å
   * @param cutoff inclusion radius R0 in Å; neighbours must be strictly closer
   * @return the sparse neighbourhood
   */
  public static StructuralNeighborhood build(float[][] distances, float cutoff) {
    if (!Float.isFinite(cutoff) || cutoff <= 0.0f) {
      throw new IllegalArgumentException("Cutoff must be finite and positive: " + cutoff);
    }
    int residues = distances.length;
    int[][] indices = new int[residues][];
    float[][] neighborDistances = new float[residues][];
    float[][] weights = new float[residues][];
    for (int residue = 0; residue < residues; residue++) {
      if (distances[residue].length != residues) {
        throw new IllegalArgumentException("Distance matrix must be square");
      }
      int count = 0;
      for (int neighbor = 0; neighbor < residues; neighbor++) {
        if (isNeighbor(residue, neighbor, distances[residue][neighbor], cutoff)) {
          count++;
        }
      }
      indices[residue] = new int[count];
      neighborDistances[residue] = new float[count];
      weights[residue] = new float[count];
      int position = 0;
      for (int neighbor = 0; neighbor < residues; neighbor++) {
        float distance = distances[residue][neighbor];
        if (isNeighbor(residue, neighbor, distance, cutoff)) {
          indices[residue][position] = neighbor;
          neighborDistances[residue][position] = distance;
          weights[residue][position] = Math.max(MIN_CONTACT_WEIGHT, 1.0f - distance / cutoff);
          position++;
        }
      }
    }
    return new StructuralNeighborhood(indices, neighborDistances, weights);
  }

  /** Returns the number of structural residues. */
  public int residueCount() {
    return residueCount;
  }

  /** Returns the structural indices of the neighbours of a residue. */
  public int[] neighborIndices(int residue) {
    return neighborIndices[residue];
  }

  /** Returns the distances to the neighbours of a residue, in the same order. */
  public float[] neighborDistances(int residue) {
    return neighborDistances[residue];
  }

  /** Returns the contact weights of the neighbours of a residue, in the same order. */
  public float[] neighborWeights(int residue) {
    return neighborWeights[residue];
  }

  private static boolean isNeighbor(int residue, int neighbor, float distance, float cutoff) {
    return neighbor != residue && Float.isFinite(distance) && distance >= 0.0f && distance < cutoff;
  }
}
