package org.transept.structure;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import org.transept.io.InvalidInputFileException;
import org.transept.io.PdbChainReader;
import org.transept.io.StructureData;
import org.transept.io.StructureDataReader;
import org.transept.model.ProteinSequence;
import org.transept.model.ProteinSet;

/**
 * Loads the structures of a protein set: distance matrices and residue numbers from
 * {@code distances/<dataset>/} and, when available, PDB files from {@code alphafold_pdb/<dataset>/}
 * to verify the residue mapping.
 *
 * <p>Missing or malformed structural files do not stop loading: the protein is reported as
 * without structure, so that callers can decide whether to accept the run.
 */
public final class StructureLoader {

  private final StructureDataReader dataReader;
  private final PdbChainReader pdbReader;
  private final SequenceStructureMapper mapper;

  /**
   * Creates a loader.
   *
   * @param policy mapping policy between sequence positions and structure rows
   */
  public StructureLoader(MappingPolicy policy) {
    this.dataReader = new StructureDataReader();
    this.pdbReader = new PdbChainReader();
    this.mapper = new SequenceStructureMapper(policy);
  }

  /**
   * Loads the structures of every protein.
   *
   * @param proteins proteins in canonical order
   * @param distanceDirectory directory with {@code <id>_D.csv} and {@code <id>_idx.csv}
   * @param pdbDirectory directory with {@code <id>.pdb}, if available
   * @return the structures and their mapping reports
   * @throws IOException if an existing file cannot be read
   */
  public ProteinStructures load(
      ProteinSet proteins, Path distanceDirectory, Optional<Path> pdbDirectory)
      throws IOException {
    List<StructureMapping> mappings = new ArrayList<>(proteins.size());
    List<float[][]> distances = new ArrayList<>(proteins.size());
    for (ProteinSequence protein : proteins.proteins()) {
      Path distanceFile = StructureDataReader.distanceFile(distanceDirectory, protein.id());
      Path indexFile = StructureDataReader.indexFile(distanceDirectory, protein.id());
      StructureMapping mapping;
      float[][] matrix = null;
      if (!Files.isRegularFile(distanceFile) || !Files.isRegularFile(indexFile)) {
        mapping =
            StructureMapping.unusable(
                protein.id(), MappingStatus.MISSING_FILES, "No files in " + distanceDirectory);
      } else {
        try {
          StructureData data = dataReader.read(distanceFile, indexFile);
          mapping = mapper.map(protein, data.residueNumbers(), pdbSequence(protein, pdbDirectory));
          matrix = data.distances();
        } catch (InvalidInputFileException exception) {
          mapping =
              StructureMapping.unusable(
                  protein.id(), MappingStatus.INVALID_FILES, exception.getMessage());
        }
      }
      mappings.add(mapping);
      distances.add(mapping.status().isUsable() ? matrix : null);
    }
    return new ProteinStructures(mappings, distances);
  }

  private Optional<String> pdbSequence(ProteinSequence protein, Optional<Path> pdbDirectory)
      throws IOException {
    Optional<String> sequence = Optional.empty();
    if (pdbDirectory.isPresent()) {
      Path pdb = pdbDirectory.get().resolve(protein.id() + ".pdb");
      if (Files.isRegularFile(pdb)) {
        sequence = Optional.of(pdbReader.read(pdb).sequence());
      }
    }
    return sequence;
  }
}
