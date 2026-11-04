# T6. Using Your Own Dataset

| | |
| --- | --- |
| **Level** | Advanced |
| **Version** | 1.0 (2026-11-04) |
| **Time** | about 30 minutes, plus predicting topologies and structures |
| **Prerequisites** | [T1. First run in five minutes](T1-first-run.md), Python 3 with NumPy and Biopython |
| **Checked with** | A copy of `classF` built with the steps below; the script of Step 3 reproduces the benchmark matrices to 0.0001 Å |

TM-M2StructAlign reads its inputs from a directory with the layout of the benchmark
([Dataset_Structural_GPCRs](https://github.com/TM-M2StructAlign/Dataset_Structural_GPCRs)).
This tutorial explains which files a dataset needs, how to produce them for your own proteins,
and how to check that the tool reads them as intended.

## Step 1: the layout

A dataset named `mydata` inside a root directory `MyRoot` needs these files:

```text
MyRoot/
├── sequences/tmregions/mydata_predicted_topologies.3line   sequences and topologies (required)
├── distances/mydata/<id>_D.csv                             distance matrix of each protein
├── distances/mydata/<id>_idx.csv                           residue numbers of its rows
├── alphafold_pdb/mydata/<id>.pdb                           structure models (recommended)
└── precomputed/mydata/<tool>_mydata.fasta                  seed alignments (at least one)
```

Then `--data-root MyRoot --dataset mydata --seeds <tool>,...` runs on it. `<id>` is the identifier
of the protein in the topology file; it must be used verbatim in the structural file names.

## Step 2: sequences and topology

The **topology file is the source of the proteins**: their identifiers, their sequences and their
order come from it. It has three lines per protein:

```text
>fzd2_human | SP+TM
MRPRSALPRLLLPLLLLPAAGPAQFHGEKGISIPDHGFCQPISIPLCTDIAYNQTIMPNLLGHTNQEDAGLEVHQFYPLV...
SSSSSSSSSSSSSSSSSSSSSSSOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO...
```

- The identifier is the text of the header before the first `|`, trimmed.
- The third line has one letter per residue. `M` marks transmembrane residues; any other letter
  (`I`, `O`, `S`, ...) is treated as non-transmembrane.

This is the 3-line format written, for example, by DeepTMHMM. A protein needs at least one `M`
residue for the transmembrane rules to matter, but the tool also accepts proteins without them.

## Step 3: structures and distance matrices

The structural objective needs, for every protein, the Cα distance matrix of a structural model
(AlphaFold2 in the benchmark) and the residue number of every row:

- `<id>_D.csv`: a square matrix of distances in Å, comma-separated, one row per line. The benchmark
  caps distances at 50 Å; only distances below the inclusion radius (15 Å) matter.
- `<id>_idx.csv`: one PDB residue number per line, in the order of the rows of the matrix.

This script computes both from the PDB models of a directory, as the benchmark does:

```python
#!/usr/bin/env python3
"""Usage: make_distances.py <pdb-directory> <output-directory>"""
import pathlib
import sys

import numpy as np
from Bio.PDB import PDBParser

MAX_DISTANCE = 50.0

source, target = map(pathlib.Path, sys.argv[1:3])
target.mkdir(parents=True, exist_ok=True)
for pdb in sorted(source.glob("*.pdb")):
    chain = next(PDBParser(QUIET=True).get_structure(pdb.stem, pdb)[0].get_chains())
    residues = [r for r in chain if r.id[0] == " " and "CA" in r]
    coordinates = np.array([r["CA"].get_coord() for r in residues], dtype=np.float64)
    distances = np.linalg.norm(coordinates[:, None, :] - coordinates[None, :, :], axis=2)
    np.savetxt(target / f"{pdb.stem}_D.csv", np.minimum(distances, MAX_DISTANCE),
               fmt="%.4f", delimiter=",")
    np.savetxt(target / f"{pdb.stem}_idx.csv", [r.id[1] for r in residues], fmt="%d")
```

```bash
python3 make_distances.py MyRoot/alphafold_pdb/mydata MyRoot/distances/mydata
```

Name the PDB files after the identifiers of the topology file (`<id>.pdb`). Keep the PDB files in
`alphafold_pdb/mydata/`: with them, the tool checks that the residues of the model match the
sequence (next step). A protein without structural files is still aligned, but does not contribute
to the lDDT.

## Step 4: how structures are matched to sequences

The residue numbers of a model do not always index the sequence: the model may be of another
isoform, or start at another residue. With the default `--mapping-policy VERIFIED`, every protein
is mapped in this order:

1. **By numbering**: residue number `n` is residue `n` of the sequence. Accepted if at least 95 %
   of the residues agree (status `NUMBERED`).
2. **By alignment**: otherwise, the sequence of the PDB model is aligned to the input sequence, and
   the aligned residues are mapped. Accepted if at least 90 % of them are identical (status
   `ALIGNED`).
3. Otherwise the structure is not used (status `LOW_IDENTITY`).

Without a PDB file, the numbering cannot be checked and is used as is (status `UNVERIFIED`).
`NUMBERED_ONLY` skips step 2, and `LEGACY` trusts the numbering without checks, as the original
implementation did.

In the benchmark, the AlphaFold models of OXER1 and MCHR1 are of a different isoform than the
GPCRdb sequences: by numbering only 8 % and 7 % of their residues agree, while the alignment maps
them with 100 % identity.

## Step 5: seed alignments

Put one or more alignments of **all** the proteins of the topology file in
`precomputed/mydata/<tool>_mydata.fasta`, in aligned FASTA. Any aligner works (MAFFT, Clustal
Omega, MUSCLE, T-Coffee, ...). The headers are matched to the identifiers of the topology file
after normalisation:

- the first word of the header is used, and lower-cased;
- UniProt headers such as `sp|P12345|NAME` give `P12345`; otherwise the first `|`-separated field
  is used;
- a chain suffix (`:A`) and a `.pdb` extension are dropped.

The ungapped sequences of a seed must be identical to the sequences of the topology file. Different
seeds give the search different starting points; identical seeds are detected and counted once.
At least **two different seeds** are needed, because the rest of the initial population is made of
their crossover children (unless the population is not larger than the number of seeds).

## Step 6: run and check

Run a short search first:

```bash
java -jar target/TM-MSAStructAligner-1.0-jar-with-dependencies.jar run \
  --data-root MyRoot --dataset mydata --seeds mafft,clustalo \
  --output-root results --run-id check --seed 1 \
  --max-evaluations 200 --population-size 20
```

Then check in `results/ejecuciones/mydata/check/metadata.json`:

- `structures`: the status of every protein. `NUMBERED` and `ALIGNED` are fine; look at the
  `identity` and `coverage` of `ALIGNED` proteins, and fix the files of any protein with
  `LOW_IDENTITY`, `MISSING_FILES` or `INVALID_FILES`. The `detail` field says what is wrong, for
  example:

  ```json
  {"protein":"smo_human","status":"MISSING_FILES","identity":0.0,"coverage":0.0,
   "detail":"No files in MyRoot/distances/mydata"}
  {"protein":"fzd1_human","status":"INVALID_FILES","identity":0.0,"coverage":0.0,
   "detail":"Distance matrix MyRoot/distances/mydata/fzd1_human_D.csv is not symmetric at (0, 5)"}
  ```

  Problems with the structural files do **not** stop the run: the protein is aligned, but the
  lDDT ignores it.
- `proteinsWithoutStructure`: the proteins that the lDDT ignores.
- `seeds`: every seed with its SHA-256, and `duplicateOf` when two seeds are identical.

Problems with the sequences or the seeds stop the run before the search: it writes `error.json`
with the reason and exits with status 1. Typical messages:

| Message (in `error.json`) | What to check |
| --- | --- |
| `Protein 'fzd2_human' not found in .../bad1_mydata.fasta` | Headers of the seed and identifiers of the topology file (Step 5) |
| `Aligned row of 'fzd2_human' does not match its sequence: residue 1 is 'W' instead of 'M'` | The seed was computed from other sequences |
| `Protein 'fzd2_human' has 565 residues but 564 topology codes` | The topology file |
| `At least two distinct seed alignments are needed to fill a population of 20` | Add another seed (Step 5) |

## Step 7: compare with a reference alignment (optional)

If a reference alignment of the dataset exists (for GPCRs, the structure-based alignments of
GPCRdb), [T4](T4-choosing-and-evaluating.md) shows how to compare the results with it.
