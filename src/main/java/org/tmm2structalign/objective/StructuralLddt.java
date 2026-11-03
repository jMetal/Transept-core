package org.tmm2structalign.objective;

import java.util.Arrays;
import org.tmm2structalign.alignment.Msa;
import org.tmm2structalign.structure.ProteinStructures;
import org.tmm2structalign.structure.StructuralNeighborhood;

/**
 * Reference-free, superposition-free structural consistency of an alignment, inspired by the
 * FoldMason MSTA lDDT; maximised.
 *
 * <p>For every ordered pair of proteins with structure (query, target) and every column where
 * both have a structurally mapped residue, the query neighbours closer than R0 are looked up in
 * the target through the alignment; each neighbour scores 1, 0.75, 0.5, 0.25 or 0 depending on
 * whether the distance difference is below 0.5, 1, 2 or 4 Å, weighted by
 * {@code 1 - d / R0}. Unaligned neighbours count in the denominator only. Per-anchor scores are
 * combined according to {@link LddtAggregation}.
 *
 * <p>Port of {@code LDDTStructuralScore}; with {@link LddtAggregation#HARMONIC_OVER_LENGTH} and
 * the same residue mapping it reproduces its values. Arithmetic is kept in {@code float} where
 * the original used it.
 */
public final class StructuralLddt implements Objective {

  /** Name used in output files. */
  public static final String NAME = "lDDT";

  /** Inclusion radius used in the published experiments, in Å. */
  public static final float DEFAULT_CUTOFF = 15.0f;

  private final ProteinStructures structures;
  private final LddtAggregation aggregation;
  private final StructuralNeighborhood[] neighborhoods;
  private final float[][][] distances;

  /**
   * Creates the objective and precomputes the neighbourhood of every structure.
   *
   * @param structures structures of the proteins, in canonical order
   * @param cutoff inclusion radius R0 in Å
   * @param aggregation how anchor scores are combined
   */
  public StructuralLddt(
      ProteinStructures structures, float cutoff, LddtAggregation aggregation) {
    this.structures = structures;
    this.aggregation = aggregation;
    this.neighborhoods = new StructuralNeighborhood[structures.size()];
    this.distances = new float[structures.size()][][];
    for (int protein = 0; protein < structures.size(); protein++) {
      distances[protein] = structures.distances(protein).orElse(null);
      neighborhoods[protein] =
          distances[protein] == null
              ? null
              : StructuralNeighborhood.build(distances[protein], cutoff);
    }
  }

  @Override
  public String name() {
    return NAME;
  }

  @Override
  public boolean isMaximized() {
    return true;
  }

  @Override
  public double evaluate(Msa msa) {
    if (msa.size() != structures.size()) {
      throw new IllegalArgumentException(
          "Alignment has " + msa.size() + " sequences but " + structures.size()
              + " structures were loaded");
    }
    double value = 0.0;
    if (msa.size() >= 2 && msa.length() > 0) {
      int length = msa.length();
      double[] columnScoreSum = new double[length];
      int[] columnScoreCount = new int[length];
      int[][] structureIndexAtColumn = structureIndexAtColumn(msa);
      int[][] columnOfStructure = columnOfStructure(structureIndexAtColumn);
      for (int query = 0; query < msa.size(); query++) {
        if (neighborhoods[query] == null) {
          continue;
        }
        for (int target = 0; target < msa.size(); target++) {
          if (query == target || distances[target] == null) {
            continue;
          }
          scoreAnchors(
              query, target, structureIndexAtColumn, columnOfStructure,
              columnScoreSum, columnScoreCount);
        }
      }
      value = aggregate(msa, columnScoreSum, columnScoreCount);
    }
    return value;
  }

  private void scoreAnchors(
      int query,
      int target,
      int[][] structureIndexAtColumn,
      int[][] columnOfStructure,
      double[] columnScoreSum,
      int[] columnScoreCount) {
    for (int column = 0; column < columnScoreSum.length; column++) {
      int queryAnchor = structureIndexAtColumn[query][column];
      int targetAnchor = structureIndexAtColumn[target][column];
      if (queryAnchor < 0 || targetAnchor < 0) {
        continue;
      }
      double score =
          directionalLddt(
              neighborhoods[query],
              distances[target],
              structureIndexAtColumn[target],
              columnOfStructure[query],
              queryAnchor,
              targetAnchor);
      if (!Double.isNaN(score)) {
        columnScoreSum[column] += score;
        columnScoreCount[column]++;
      }
    }
  }

  private static double directionalLddt(
      StructuralNeighborhood queryNeighborhood,
      float[][] targetDistances,
      int[] targetStructureAtColumn,
      int[] queryColumnOfStructure,
      int queryAnchor,
      int targetAnchor) {
    double result = Double.NaN;
    if (queryAnchor < queryNeighborhood.residueCount()
        && targetAnchor < targetDistances.length
        && queryNeighborhood.neighborIndices(queryAnchor).length > 0) {
      int[] neighborIndices = queryNeighborhood.neighborIndices(queryAnchor);
      float[] neighborDistances = queryNeighborhood.neighborDistances(queryAnchor);
      float[] neighborWeights = queryNeighborhood.neighborWeights(queryAnchor);
      double weightedScoreSum = 0.0;
      double weightSum = 0.0;
      for (int position = 0; position < neighborIndices.length; position++) {
        float queryWeight = neighborWeights[position];
        weightSum += queryWeight;
        int column = queryColumnOfStructure[neighborIndices[position]];
        if (column < 0) {
          continue;
        }
        int targetNeighbor = targetStructureAtColumn[column];
        if (targetNeighbor < 0 || targetNeighbor >= targetDistances[targetAnchor].length) {
          continue;
        }
        float targetDistance = targetDistances[targetAnchor][targetNeighbor];
        if (Float.isFinite(targetDistance)) {
          float difference = Math.abs(neighborDistances[position] - targetDistance);
          weightedScoreSum += queryWeight * thresholdScore(difference);
        }
      }
      result = weightSum > 0.0 ? weightedScoreSum / weightSum : Double.NaN;
    }
    return result;
  }

  private static double thresholdScore(float difference) {
    double score;
    if (difference >= 4.0f) {
      score = 0.0;
    } else if (difference >= 2.0f) {
      score = 0.25;
    } else if (difference >= 1.0f) {
      score = 0.50;
    } else if (difference >= 0.5f) {
      score = 0.75;
    } else {
      score = 1.0;
    }
    return score;
  }

  private double aggregate(Msa msa, double[] columnScoreSum, int[] columnScoreCount) {
    double value;
    if (aggregation == LddtAggregation.PAIR_MEAN) {
      double sum = Arrays.stream(columnScoreSum).sum();
      long count = Arrays.stream(columnScoreCount).asLongStream().sum();
      value = count == 0 ? 0.0 : sum / count;
    } else {
      double total = 0.0;
      for (int column = 0; column < msa.length(); column++) {
        double occupancy = msa.occupancy(column) / (double) msa.size();
        double structuralScore =
            columnScoreCount[column] == 0
                ? 0.0
                : columnScoreSum[column] / columnScoreCount[column];
        if (occupancy > 0.0 && structuralScore > 0.0) {
          total += 2.0 * structuralScore * occupancy / (structuralScore + occupancy);
        }
      }
      value = total / msa.length();
    }
    return value;
  }

  /** Structural row at every column, {@code -1} for gaps and unmapped residues. */
  private int[][] structureIndexAtColumn(Msa msa) {
    int[][] result = new int[msa.size()][msa.length()];
    for (int protein = 0; protein < msa.size(); protein++) {
      for (int column = 0; column < msa.length(); column++) {
        int residue = msa.residueIndex(protein, column);
        result[protein][column] =
            residue < 0 || distances[protein] == null
                ? -1
                : structures.structureIndex(protein, residue);
      }
    }
    return result;
  }

  /** First column holding every structural row, {@code -1} if the row is not aligned. */
  private int[][] columnOfStructure(int[][] structureIndexAtColumn) {
    int[][] result = new int[structureIndexAtColumn.length][];
    for (int protein = 0; protein < structureIndexAtColumn.length; protein++) {
      int rows = distances[protein] == null ? 0 : distances[protein].length;
      result[protein] = new int[rows];
      Arrays.fill(result[protein], -1);
      for (int column = 0; column < structureIndexAtColumn[protein].length; column++) {
        int row = structureIndexAtColumn[protein][column];
        if (row >= 0 && row < rows && result[protein][row] < 0) {
          result[protein][row] = column;
        }
      }
    }
    return result;
  }
}
